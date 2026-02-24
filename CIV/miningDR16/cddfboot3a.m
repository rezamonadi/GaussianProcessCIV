%% C IV CDDF + β Slope + Ω_C IV  — Final Clean Version
clearvars; close all; clc;

%% -------------------------------------------------------------------------
% 1. LOAD DATA
% -------------------------------------------------------------------------
load sampled_logN_per_absorber.mat          % From posterior sampling
load processed_qsos_dr16_parameters.mat     % Includes savingCat structure

%% -------------------------------------------------------------------------
% 2. Collapse Each Absorber Posterior to a Single log N (mode and mean)
% -------------------------------------------------------------------------
edges_mode = 12.5:0.05:16.1;
n_abs = numel(sampledAbsorbers);
logN_mean = zeros(n_abs, 1);
logN_mode = zeros(n_abs, 1);

for i = 1:n_abs
    samples = sampledAbsorbers(i).logN_samples;
    logN_mean(i) = mean(samples);
    logN_mode(i) = mode(samples);
end

% Filter out absorbers with logN at sampling grid edges (unconstrained)
is_edge = logN_mode <= 12.50 | logN_mode >= 16.1;
fprintf('⚠️  Excluding %d absorbers at logN grid edges (12.5 or 16.1)\n', sum(is_edge));
logN_mode = logN_mode(~is_edge);
logN_mean = logN_mean(~is_edge);
sampledAbsorbers = sampledAbsorbers(~is_edge);

%% -------------------------------------------------------------------------
% 3. Define CDDF Bins and Completeness Correction
% -------------------------------------------------------------------------
logN_edges = 11:0.25:19;
dlogN = diff(logN_edges(1:2));
logN_centers = logN_edges(1:end-1) + dlogN/2;
N_centers = 10.^logN_centers;
n_bins = numel(logN_centers);

completion = @(logN) min(max((logN - 13.0) / (14.2 - 13.0), 0.3), 1.0);
completeness_correction = completion(logN_centers);

%% -------------------------------------------------------------------------
% 4. Compute Absorption Path ΔX
% -------------------------------------------------------------------------
Omega_M = 0.3; Omega_Lambda = 0.7;
H0 = 70; c = 2.99792458e5;                     % km/s
H0_cgs = H0 * 1e5 / 3.085677581e24;
c_cgs = c * 1e5;                               % cm/s

dX_dz = @(z) (1+z).^2 ./ sqrt(Omega_M*(1+z).^3 + Omega_Lambda);

used_qsos = unique([sampledAbsorbers.qso_idx]);
z_min_all = savingCat.all_min_z_c4s(used_qsos);
z_max_all = savingCat.all_max_z_c4s(used_qsos);

lambda_min_obs = 3600; lambda_max_obs = 10400; % SDSS spectral range
lambda_rest_c4 = 1550;
z_c4_min_allowed = lambda_min_obs / lambda_rest_c4 - 1;
z_c4_max_allowed = lambda_max_obs / lambda_rest_c4 - 1;

z_min_clipped = max(z_min_all, z_c4_min_allowed);
z_max_clipped = min(z_max_all, z_c4_max_allowed);

dz = z_max_clipped - z_min_clipped;
min_dz_threshold = 0;
valid_mask = dz >= min_dz_threshold & isfinite(z_min_clipped) & isfinite(z_max_clipped);

z_min_valid = z_min_clipped(valid_mask);
z_max_valid = z_max_clipped(valid_mask);

n_valid = numel(z_min_valid);
delta_X_sightlines = arrayfun(@(a,b) trapz(linspace(a,b,200), dX_dz(linspace(a,b,200))), z_min_valid, z_max_valid);
total_dX = nansum(delta_X_sightlines);
fprintf('Using %d sightlines | Total ΔX = %.2f\n', n_valid, total_dX);

%% -------------------------------------------------------------------------
% 5. Bootstrap CDDF and Fit β
% -------------------------------------------------------------------------
n_trials = 200;
f_trials = zeros(n_trials, n_bins);
beta_vals = zeros(n_trials,1);

n_abs = numel(sampledAbsorbers);
pre_sampled_logN = arrayfun(@(s) s.logN_samples(randi(numel(s.logN_samples))), sampledAbsorbers);

valid_mask = isfinite(pre_sampled_logN);
pre_sampled_logN = pre_sampled_logN(valid_mask);

for t = 1:n_trials
    idx = randsample(length(pre_sampled_logN), length(pre_sampled_logN), true);
    logN_trial = pre_sampled_logN(idx);
    counts = histcounts(logN_trial, logN_edges);
    f = counts ./ (dlogN * total_dX);
    f = f ./ (N_centers .* log(10));
    f = f ./ completeness_correction;
    f_trials(t,:) = f;

    logf_trial = log10(f);
    fit_mask = logN_centers >= 13.4 & logN_centers <= 14.8 & isfinite(logf_trial);
    if sum(fit_mask) >= 2
        p = polyfit(logN_centers(fit_mask), logf_trial(fit_mask), 1);
        beta_vals(t) = -p(1);
    else
        beta_vals(t) = NaN;
    end
end

f_mean = mean(f_trials,1);
logf = log10(f_mean);
beta = nanmean(beta_vals);
beta_std = nanstd(beta_vals);
fprintf('Best-fit β = %.3f ± %.3f\n', beta, beta_std);

%% -------------------------------------------------------------------------
% 6. Plot CDDF + β Fit with Literature Comparison
% -------------------------------------------------------------------------
figure;
errorbar(logN_centers, logf, log10(f_mean)-log10(min(f_trials)), log10(max(f_trials))-log10(f_mean), ...
         'k^','MarkerFaceColor','k','LineStyle','none','CapSize',0);
xlabel('$\log_{10}(N_{\mathrm{C\,IV}}/\mathrm{cm}^{-2})$', 'Interpreter','latex');
ylabel('$\log_{10}[f(N,X)]$', 'Interpreter','latex');
title('C IV Column-Density Distribution Function'); grid on;

hold on;
fit_mask = logN_centers >= 13.4 & logN_centers <= 14.8 & isfinite(logf);
p = polyfit(logN_centers(fit_mask), logf(fit_mask), 1);
y_line = polyval(p, logN_centers(fit_mask));
plot(logN_centers(fit_mask), y_line, 'r--', 'LineWidth', 1.5);
text(mean(logN_centers(fit_mask)), mean(logf(fit_mask)), ...
     string(sprintf('$\\beta = %.3f$', beta)), ...
     'Color','r', 'FontSize', 12, 'Interpreter','latex', ...
     'HorizontalAlignment','center');

logN_song = [12.5, 13.0, 13.5, 14.0, 14.5];
logf_song = [-11.5, -12.2, -13.0, -13.8, -14.6];
err_song =  [0.2,    0.2,   0.2,   0.3,    0.4];

logN_cooksey = [13.3, 13.6, 13.9, 14.2];
logf_cooksey = [-13.0, -13.4, -13.9, -14.5];
err_cooksey  = [0.1,    0.1,   0.15,  0.2];

logN_dodorico = [13.0, 13.5, 14.0, 14.5];
logf_dodorico = [-12.3, -13.0, -13.8, -14.7];
err_dodorico  = [0.2,    0.2,   0.2,   0.25];

errorbar(logN_song, logf_song, err_song, 'mo', 'LineWidth', 1.2, 'DisplayName','Songaila 2001');
errorbar(logN_cooksey, logf_cooksey, err_cooksey, 'bs', 'LineWidth', 1.2, 'DisplayName','Cooksey et al. 2013');
errorbar(logN_dodorico, logf_dodorico, err_dodorico, 'gd', 'LineWidth', 1.2, 'DisplayName','D''Odorico et al. 2010');
legend('Location','southwest');
%% -------------------------------------------------------------------------
% 7. Ω_C IV from logN = 12.5.0–15.0 range
% -------------------------------------------------------------------------

logN_lo = 12.5;      % or 12.0
logN_hi = 15.0;
valid_abs_mask = logN_mean >= logN_lo & logN_mean <= logN_hi;
logN_fit = logN_mean(valid_abs_mask);
comp_fit = completion(logN_fit);


% Physical constants
A_C   = 12.0107;
m_u   = 1.66053906660e-24;
m_CIV = A_C * m_u;
G     = 6.67430e-8;
rho_crit = 3 * H0_cgs^2 / (8 * pi * G);
prefac = (H0_cgs * m_CIV) / (c_cgs * rho_crit);

% Compute Omega_C IV
N_vals_corr = 10.^logN_fit ./ comp_fit;
Omega_CIV_corr = prefac * sum(N_vals_corr) / total_dX;
fprintf('Ω_C IV (corrected, %.1f ≤ log N ≤ %.1f) = %.3e\n', logN_lo, logN_hi, Omega_CIV_corr);

% Bootstrap Omega_C IV
nboot = 500;
Omega_boot = zeros(nboot, 1);
for b = 1:nboot
    idx = randsample(numel(logN_fit), numel(logN_fit), true);
    N_b = 10.^logN_fit(idx) ./ comp_fit(idx);
    Omega_boot(b) = prefac * sum(N_b) / total_dX;
end

Omega_mean = mean(Omega_boot);
Omega_std  = std(Omega_boot);
fprintf('Ω_C IV = %.3e ± %.3e (bootstrap 1σ, %d trials)\n', Omega_mean, Omega_std, nboot);
%% -------------------------------------------------------------------------
% 8. Redshift-Binned CDDF Slope (β), Ω_C IV, and dN/dX
% -------------------------------------------------------------------------

z_edges = [1.5, 2.0, 2.5, 3.0, 4.0, 5.5];
z_centers = (z_edges(1:end-1) + z_edges(2:end)) / 2;
nbins = numel(z_centers);

beta_z = nan(nbins,1);
omega_z = nan(nbins,1);
omega_err = nan(nbins,1);
dN_dX_z = nan(nbins,1);

fprintf('\n=== Redshift-Binned CDDF Summary ===\n');

logN_mode_all = arrayfun(@(s) mode(s.logN_samples), sampledAbsorbers)';
z_abs_all = savingCat.all_map_z_c4L2([sampledAbsorbers.qso_idx]);
z_min_all = savingCat.all_min_z_c4s(used_qsos);
z_max_all = savingCat.all_max_z_c4s(used_qsos);

for b = 1:nbins
    zmin = z_edges(b);
    zmax = z_edges(b+1);

    in_bin = z_abs_all >= zmin & z_abs_all < zmax;
    logN_bin = logN_mode_all(in_bin);
    if isempty(logN_bin), continue; end

    contributing = z_max_all > zmin & z_min_all < zmax;
    zmin_clip = max(z_min_all(contributing), zmin);
    zmax_clip = min(z_max_all(contributing), zmax);

    delta_X = arrayfun(@(a,b) trapz(linspace(a,b,200), dX_dz(linspace(a,b,200))), ...
                       zmin_clip, zmax_clip);
    total_dX_bin = nansum(delta_X);

    % Histogram and CDDF
    counts = histcounts(logN_bin, logN_edges);
    f = counts ./ (dlogN * total_dX_bin);
    f = f ./ (N_centers .* log(10));
    f = f ./ completeness_correction;

    fit_mask = logN_centers >= 13.4 & logN_centers <= 14.8 & isfinite(f) & f > 0;
    if sum(fit_mask) >= 2
        p = polyfit(logN_centers(fit_mask), log10(f(fit_mask)), 1);
        beta_z(b) = -p(1);
    end

    % Ω_C IV
    omega_mask = logN_bin >= logN_lo & logN_bin <= logN_hi;
    logN_fit_bin = logN_bin(omega_mask);
    comp_bin = completion(logN_fit_bin);
    N_vals_corr = 10.^logN_fit_bin ./ comp_bin;
    omega_z(b) = prefac * sum(N_vals_corr) / total_dX_bin;

    % Bootstrap Ω_C IV
    nboot = 200;
    omega_boot = zeros(nboot,1);
    for i = 1:nboot
        idx = randsample(numel(logN_fit_bin), numel(logN_fit_bin), true);
        N_b = 10.^logN_fit_bin(idx) ./ comp_bin(idx);
        omega_boot(i) = prefac * sum(N_b) / total_dX_bin;
    end
    omega_err(b) = std(omega_boot);

    % dN/dX
    dN_dX_z(b) = sum(omega_mask) / total_dX_bin;

    fprintf('z = %.1f–%.1f: β = %.2f, Ω_C IV = %.2e ± %.2e, dN/dX = %.3f\n', ...
            zmin, zmax, beta_z(b), omega_z(b), omega_err(b), dN_dX_z(b));
end

% Plot β(z)
figure;
plot(z_centers, beta_z, 'r-o', 'LineWidth', 1.5);
xlabel('Redshift'); ylabel('\beta');
title('CDDF Slope \beta vs. Redshift');
grid on;

% Plot Ω_C IV(z)
figure;
errorbar(z_centers, omega_z, omega_err, 'ko-', 'MarkerFaceColor','k');
xlabel('Redshift'); ylabel('$\Omega_{\mathrm{C\,IV}}$', 'Interpreter','latex');
title('$\Omega_{\mathrm{C\,IV}}$ vs. Redshift', 'Interpreter','latex');
grid on;

% Plot dN/dX(z)
figure;
bar(z_centers, dN_dX_z, 'FaceColor', [0.5 0.5 1]);
xlabel('Redshift'); ylabel('dN/dX');
title('C IV Absorber Incidence vs. Redshift');
grid on;
%% -------------------------------------------------------------------------
% 9. Plot CDDFs by Redshift Bin
% -------------------------------------------------------------------------
colors = lines(length(z_edges)-1);  % Assign distinct colors to each redshift bin
figure; hold on;
% Precompute ΔX per redshift bin (X_per_bin)
X_per_bin = zeros(length(z_edges)-1, 1);

for b = 1:length(z_edges)-1
    zmin = z_edges(b);
    zmax = z_edges(b+1);

    % Identify contributing sightlines
    contributing = z_max_all > zmin & z_min_all < zmax;
    zmin_clip = max(z_min_all(contributing), zmin);
    zmax_clip = min(z_max_all(contributing), zmax);

    % Compute ΔX per sightline and sum
    dX_vals = arrayfun(@(a,b) trapz(linspace(a,b,200), dX_dz(linspace(a,b,200))), ...
                       zmin_clip, zmax_clip);
    X_per_bin(b) = nansum(dX_vals);
end

for b = 1:length(z_edges)-1
    zmin = z_edges(b);
    zmax = z_edges(b+1);

    in_bin = z_abs_all >= zmin & z_abs_all < zmax;
    logN_bin = logN_mode_all(in_bin);
    if numel(logN_bin) < 50, continue; end  % Skip sparse bins

    counts = histcounts(logN_bin, logN_edges);
    f = counts ./ (dlogN * X_per_bin(b));
    f = f ./ (N_centers .* log(10));
    f = f ./ completeness_correction;

    plot(logN_centers, log10(f), '-', 'Color', colors(b,:), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('z = %.1f–%.1f', zmin, zmax));
end
plot(logN_centers, log10(f_mean), 'k-', 'LineWidth', 2, 'DisplayName', 'Global CDDF');

xlabel('$\log_{10}(N_{\mathrm{C\,IV}} / \mathrm{cm}^{-2})$', 'Interpreter','latex');
ylabel('$\log_{10}[f(N,X)]$', 'Interpreter','latex');
title('Redshift-Binned C IV CDDF');
legend('Location','southwest');
grid on;
% -------------------------------------------------------------------------
% 10. Fit Omega_CIV(z) = A * (1 + z)^gamma
% -------------------------------------------------------------------------

% Filter to bins with finite Ω values
valid = isfinite(omega_z) & omega_z > 0;

% Fit in log-log space
z_fit = z_centers(valid);
log1pz = log10(1 + z_fit);
logOmega = log10(omega_z(valid));
logOmega_err = omega_err(valid) ./ (omega_z(valid) * log(10));  % error in log10(Ω)

% Weighted linear fit: log10(Ω) = log10(A) + γ * log10(1+z)
p = polyfit(log1pz, logOmega, 1);
gamma = p(1);
logA = p(2);
A = 10^logA;

% Bootstrap fit uncertainty (optional)
nboot = 500;
gamma_boot = zeros(nboot,1);
for b = 1:nboot
    idx = randsample(numel(log1pz), numel(log1pz), true);
    p_b = polyfit(log1pz(idx), logOmega(idx), 1);
    gamma_boot(b) = p_b(1);
end
gamma_err = std(gamma_boot);

% Report to console (no LaTeX)
fprintf('Omega_CIV(z) = A * (1 + z)^gamma\n');
fprintf('Best-fit A = %.2e, gamma = %.2f ± %.2f\n', A, gamma, gamma_err);

% Plot fit
z_plot = linspace(min(z_fit), max(z_fit), 200);
Omega_fit = A * (1 + z_plot).^gamma;

figure;
errorbar(z_centers, omega_z, omega_err, 'ko', 'MarkerFaceColor','k'); hold on;
plot(z_plot, Omega_fit, 'r-', 'LineWidth', 2);
xlabel('Redshift');
ylabel('$\Omega_{\mathrm{C\,IV}}$', 'Interpreter','latex');
title('$\Omega_{\mathrm{C\,IV}}(z)$ fit: $A(1+z)^\gamma$', 'Interpreter','latex');
legend("Binned Data", ...
       sprintf("$\\gamma = %.2f \\pm %.2f$", gamma, gamma_err), ...
       'Interpreter','latex', 'Location','northeast');
grid on;
% -------------------------------------------------------------------------
% 11. Overlay literature & simulation Ω_CIV(z) on fitted evolution plot
% -------------------------------------------------------------------------

% Simulated and literature Ω_CIV(z) curves (linear units)
z_lit = [1.75, 2.25, 2.75, 3.25, 4.00, 4.75];

% Approximate values from publications
omega_dodorico = [9e-9, 7e-9, 5.5e-9, 4e-9, 2.5e-9, 2e-9];    % D'Odorico+2010
omega_simcoe   = [1.5e-8, 1.2e-8, 8e-9, 6e-9, 4e-9, 3e-9];     % Simcoe+2011
omega_fire     = [1.2e-8, 9e-9, 6e-9, 4e-9, 2.5e-9, 1.5e-9];   % FIRE (Hafen+2019)

% Create the fitted curve again if needed
z_plot = linspace(min(z_fit), max(z_fit), 200);
Omega_fit = A * (1 + z_plot).^gamma;

% Plot all on the same figure
figure;
errorbar(z_centers, omega_z, omega_err, 'ko', 'MarkerFaceColor','k'); hold on;
plot(z_plot, Omega_fit, 'r-', 'LineWidth', 2, 'DisplayName', 'Best-fit');

% Overlay literature and simulation tracks
plot(z_lit, omega_dodorico, 'g--o', 'LineWidth', 1.2, 'DisplayName', 'D''Odorico+2010');
plot(z_lit, omega_simcoe,   'b--s', 'LineWidth', 1.2, 'DisplayName', 'Simcoe+2011');
plot(z_lit, omega_fire,     'm--^', 'LineWidth', 1.2, 'DisplayName', 'FIRE (Hafen+2019)');

% Label and legend
xlabel('Redshift');
ylabel('$\Omega_{\mathrm{C\,IV}}$', 'Interpreter','latex');
title('$\Omega_{\mathrm{C\,IV}}(z)$ vs. Literature and Simulations', 'Interpreter','latex');
legend('Interpreter','latex', 'Location','northeast');
grid on;

%% -------------------------------------------------------------------------
% 6. Plot CDDF + β Fit using logN_mean Bins (Replacing logN_centers)
% -------------------------------------------------------------------------

% Bin logN_mean into same edges
counts = histcounts(logN_mean, logN_edges);
logN_mean_bins = zeros(size(logN_centers));
for i = 1:n_bins
    in_bin = logN_mean >= logN_edges(i) & logN_mean < logN_edges(i+1);
    if any(in_bin)
        logN_mean_bins(i) = mean(logN_mean(in_bin));
    else
        logN_mean_bins(i) = NaN;
    end
end

% Compute CDDF from f_mean (already calculated from bootstraps in Section 5)
logf = log10(f_mean);

% Only keep bins with valid centers and logf
valid_idx = ~isnan(logN_mean_bins) & isfinite(logf);
logN_mean_valid = logN_mean_bins(valid_idx);
logf_valid = logf(valid_idx);

% Plot CDDF
figure;
errorbar(logN_mean_valid, logf_valid, ...
         log10(f_mean(valid_idx)) - log10(min(f_trials(:,valid_idx),[],1)), ...
         log10(max(f_trials(:,valid_idx),[],1)) - log10(f_mean(valid_idx)), ...
         'k^','MarkerFaceColor','k','LineStyle','none','CapSize',0);
xlabel('$\log_{10}(N_{\mathrm{C\,IV}}/\mathrm{cm}^{-2})$', 'Interpreter','latex');
ylabel('$\log_{10}[f(N,X)]$', 'Interpreter','latex');
title('C IV Column-Density Distribution Function (Binned by logN\_mean)');
grid on; hold on;

% Fit β in range 13.4 ≤ logN ≤ 14.8
fit_mask = logN_mean_valid >= 13.4 & logN_mean_valid <= 14.8;
p = polyfit(logN_mean_valid(fit_mask), logf_valid(fit_mask), 1);
y_line = polyval(p, logN_mean_valid(fit_mask));
plot(logN_mean_valid(fit_mask), y_line, 'r--', 'LineWidth', 1.5);
text(mean(logN_mean_valid(fit_mask)), mean(logf_valid(fit_mask)), ...
     string(sprintf('$\\beta = %.3f$', beta)), ...
     'Color','r', 'FontSize', 12, 'Interpreter','latex', ...
     'HorizontalAlignment','center');

% Literature comparison
logN_song = [12.5, 13.0, 13.5, 14.0, 14.5];
logf_song = [-11.5, -12.2, -13.0, -13.8, -14.6];
err_song =  [0.2,    0.2,   0.2,   0.3,    0.4];

logN_cooksey = [13.3, 13.6, 13.9, 14.2];
logf_cooksey = [-13.0, -13.4, -13.9, -14.5];
err_cooksey  = [0.1,    0.1,   0.15,  0.2];

logN_dodorico = [13.0, 13.5, 14.0, 14.5];
logf_dodorico = [-12.3, -13.0, -13.8, -14.7];
err_dodorico  = [0.2,    0.2,   0.2,   0.25];

errorbar(logN_song, logf_song, err_song, 'mo', 'LineWidth', 1.2, 'DisplayName','Songaila 2001');
errorbar(logN_cooksey, logf_cooksey, err_cooksey, 'bs', 'LineWidth', 1.2, 'DisplayName','Cooksey et al. 2013');
errorbar(logN_dodorico, logf_dodorico, err_dodorico, 'gd', 'LineWidth', 1.2, 'DisplayName','D''Odorico et al. 2010');
legend('Location','southwest');
