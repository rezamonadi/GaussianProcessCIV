%% C IV Likelihood Sampling Grid Extraction (full column density range)
clearvars; close all; 

load('processed_qsos_dr16_parameters.mat');
m = matfile('processed_qsos_dr16_Likelihood.mat');

%% Assume HDF5 shape: [QSO, 10000 weights, 7 absorbers]
info = h5info(m.Properties.Source, '/LikelihoodCat/all_sample_log_likelihoods_c4L2');
[~, nweights, ~] = deal(info.Dataspace.Size(1), info.Dataspace.Size(2), info.Dataspace.Size(3));

%% Define the column density grid range explicitly
logN_grid = linspace(12.5, 16.5, nweights);  % Full grid: 10000 values from logN = 12 to 17

%% Save the grid to .mat file
save('logN_grid.mat', 'logN_grid');

fprintf('✅ Saved logN_grid with %d values spanning [%.1f, %.1f] to logN_grid.mat\n', ...
        numel(logN_grid), logN_grid(1), logN_grid(end));

% %%––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––%%
% %  C IV CDDF with bootstrap errors using all 7 absorbers per QSO   %
% %%––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––%%
% 
% clearvars; close all; clc;
% load('processed_qsos_dr16_parameters.mat');
% m = matfile('processed_qsos_dr16_Likelihood.mat');
% 
% %% Cosmology & ΔX calculation
% Omega_M = 0.3; Omega_Lambda = 0.7; c = 299792.458; H0 = 70;
% dX_dz = @(z) (1+z).^2 ./ sqrt(Omega_M*(1+z).^3 + Omega_Lambda);
% min_z = savingCat.all_min_z_c4s(:); max_z = savingCat.all_max_z_c4s(:);
% n_qsos = numel(min_z); delta_X = zeros(n_qsos,1);
% for i = 1:n_qsos
%     if min_z(i) < max_z(i)
%         zgrid = linspace(min_z(i), max_z(i), 200);
%         delta_X(i) = trapz(zgrid, dX_dz(zgrid));
%     end
% end
% total_dX = sum(delta_X);
% 
% %% Bin definitions
% logN_edges   = 12:0.1:17;
% dlogN        = diff(logN_edges(1:2));
% logN_centers = logN_edges(1:end-1) + dlogN/2;
% n_bins       = numel(logN_centers);
% N_centers    = 10.^logN_centers;
% 
% %% Bootstrap setup
% nboot       = 200;
% f_N_X_boot  = zeros(nboot, n_bins);
% 
% %% HDF5 shape: [QSO, 10000 weights, 7 absorbers]
% info = h5info(m.Properties.Source, '/LikelihoodCat/all_sample_log_likelihoods_c4L2');
% [nqso_h5, nweights, max_abs_per_qso] = deal(info.Dataspace.Size(1), info.Dataspace.Size(2), info.Dataspace.Size(3));
% 
% %% Access maps
% logN_matrix = savingCat.all_map_N_c4L2;
% p_matrix    = savingCat.all_p_c4;
% 
% %% Loop over all QSOs
% for qso_idx = 1:n_qsos
%     q_start = tic;
% 
%     logN_row = logN_matrix(qso_idx, :);
%     p_row    = p_matrix(qso_idx, :);
%     valid_absorbers = find(~isnan(logN_row) & logN_row > 0 & p_row > 0.95);
% 
%     if isempty(valid_absorbers)
%         fprintf('Skipped QSO %d (no valid absorbers)\n', qso_idx);
%         continue;
%     end
% 
%     % Loop over all valid absorbers
%     for a = 1:length(valid_absorbers)
%         abs_local_idx = valid_absorbers(a);  % 1–7
%         logN_val = logN_row(abs_local_idx);
% 
%         if qso_idx > nqso_h5 || abs_local_idx > max_abs_per_qso
%             fprintf('QSO %d absorber %d: out of bounds, skipped\n', qso_idx, abs_local_idx);
%             continue;
%         end
% 
%         try
%             weights = h5read(m.Properties.Source, ...
%                 '/LikelihoodCat/all_sample_log_likelihoods_c4L2', ...
%                 [qso_idx, 1, abs_local_idx], [1, nweights, 1]);
% 
%             if isempty(weights)
%                 fprintf('QSO %d absorber %d: empty weights, skipped\n', qso_idx, abs_local_idx);
%                 continue;
%             end
%         catch ME
%             fprintf('QSO %d absorber %d: HDF5 read failed (%s)\n', qso_idx, abs_local_idx, ME.message);
%             continue;
%         end
% 
%         weights = reshape(weights, [], 1);
%         if any(isnan(weights)) || all(weights <= 0)
%             fprintf('QSO %d absorber %d: invalid weights, skipped\n', qso_idx, abs_local_idx);
%             continue;
%         end
% 
%         weights = weights / sum(weights);
% 
%         for b = 1:nboot
%             sampled_logN = logN_val;  % Replace with weighted draw if needed
%             bin_idx = find(logN_edges <= sampled_logN, 1, 'last');
%             if bin_idx < n_bins
%                 f_N_X_boot(b, bin_idx) = f_N_X_boot(b, bin_idx) + 1;
%             end
%         end
%     end
% 
%     fprintf('Done QSO %d/%d in %.2f sec\n', qso_idx, n_qsos, toc(q_start));
%     drawnow; pause(0.01);
% 
%     % Checkpoint save every 10000 QSOs
%     if mod(qso_idx, 10000) == 0
%         save(sprintf('checkpoint_cddf_qso%d.mat', qso_idx), ...
%             'f_N_X_boot', 'logN_edges', 'N_centers', 'total_dX', 'nboot');
%         fprintf('Checkpoint saved at QSO %d\n', qso_idx);
%     end
% end
% 
% %% Normalize CDDF
% f_N_X_boot = f_N_X_boot ./ (dlogN * total_dX);
% f_N_X_boot = f_N_X_boot ./ (N_centers .* log(10));
% f_N_X_mean = mean(f_N_X_boot, 1);
% f_N_X_std  = std(f_N_X_boot, 0, 1);
% logf_N_X_mean = log10(f_N_X_mean);
% sigma_logf    = f_N_X_std ./ (f_N_X_mean * log(10));
% 
% %% Plot CDDF
% figure('Position',[100 100 600 450]);
% errorbar(logN_centers, logf_N_X_mean, sigma_logf, ...
%     'Marker', '^', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', ...
%     'MarkerSize', 2, 'LineStyle', 'none', 'CapSize', 1);
% xlabel('$\\log_{10}\\left[N_{\\mathrm{C\\,IV}} / \\mathrm{cm}^{-2}\\right]$', ...
%     'Interpreter', 'latex', 'FontSize', 12);
% ylabel('$\\log_{10}\\left[f(N,X)\\right]$', ...
%     'Interpreter', 'latex', 'FontSize', 12);
% grid on;
% set(gca, 'GridLineStyle', '-', 'GridAlpha', 0.2);
