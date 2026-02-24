% --- setup ------------------------------------------------------
% original catalog
load('processed_qsos_dr16_parameters.mat');  
n_qsos = numel(savingCat.all_min_z_c4s);

% open the likelihood file
mL    = matfile('processed_qsos_dr16_Likelihood.mat');
src   = mL.Properties.Source;
nSamp = 10000;          % number of posterior samples
Mdraw = 100;            % weighted draws per absorber
pthr  = 0.85;           % or 0.95 if you prefer

% --- find grouping of absorbers into quasars -------------------
% assume savingCat gives you the number of c4 absorbers per quasar:
% (if your field is named differently, swap in the right one)
nAbsPerQSO = savingCat.all_map_N_c4L2(:);

% compute the global indices of each absorber
startIdx = [1; cumsum(nAbsPerQSO(1:end-1))+1];
endIdx   = cumsum(nAbsPerQSO);

% --- loop over every quasar & its absorbers -------------------
logN_weighted_all = [];  % will collect *all* draws

for iq = 1:n_qsos
  for j = 1:nAbsPerQSO(iq)
    k_global = startIdx(iq) + (j-1);      % index into your 1D arrays
    if savingCat.all_p_c4(k_global) < pthr
      continue   % skip low-prob absorbers
    end

    % 1) load log-likelihood samples for this (iq,j)
    %    dataset dims are [quasar, samples, absorber]
    logL = h5read(src, ...
      '/LikelihoodCat/all_sample_log_likelihoods_c4L2', ...
      [iq, 1, j], [1, nSamp, 1]);
    logL = squeeze(logL);

    % 2) turn into normalized weights
    w = exp(logL - max(logL));    % avoid overflow
    w = w / sum(w);

    % 3) draw M samples *with* those weights
    idx = randsample(nSamp, Mdraw, true, w);

    % 4) load the corresponding column-density samples
    %    (swap the dataset name if yours is different!)
    Nsam = h5read(src, ...
      '/LikelihoodCat/all_sample_map_N_c4L2', ...
      [iq, 1, j], [1, nSamp, 1]);
    Nsam = squeeze(Nsam);

    % 5) collect the weighted log10(N) draws
    logN_weighted_all = [ logN_weighted_all; log10(Nsam(idx)) ];
  end
end

% --- rebuild your histogram & CDDF with these draws -------------
% use exactly the same edges & ΔlogN as before
logN_edges   = 12:0.025:17;
logN_centers = logN_edges(1:end-1) + 0.0125;
dlogN        = logN_edges(2)-logN_edges(1);
counts_wt    = histcounts(logN_weighted_all, logN_edges);

% total ΔX from before (you’ve already computed total_dX)
f_logN_X_wt  = counts_wt / (dlogN * total_dX);
N_centers    = 10.^logN_centers;
f_N_X_wt     = f_logN_X_wt ./ (N_centers * log(10));

% plot
figure;
plot(logN_centers, log10(f_N_X_wt), 'ro');
xlabel('log_{10}(N_{C IV})');   ylabel('log_{10}[f(N,X)]');
title(sprintf('CDDF with p>%.2f weights (M=%d draws)',pthr,Mdraw));
grid on;
