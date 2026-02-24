%% C IV Absorber Sampling Script — Final Fixed Version (log-likelihoods)
clearvars; close all; 

load('processed_qsos_dr16_parameters.mat');
load('logN_grid.mat');  % Must be 1×10000 double
fprintf('logN_grid size: %s\n', mat2str(size(logN_grid)));

m = matfile('processed_qsos_dr16_Likelihood.mat');
info = h5info(m.Properties.Source, '/LikelihoodCat/all_sample_log_likelihoods_c4L2');
[nqso_h5, nweights, max_abs_per_qso] = deal(info.Dataspace.Size(1), info.Dataspace.Size(2), info.Dataspace.Size(3));

logN_matrix = savingCat.all_map_N_c4L2;
p_matrix    = savingCat.all_p_c4;

sampledAbsorbers = struct('qso_idx', {}, 'abs_idx', {}, 'logN_samples', {});
reasonsSkipped = struct('qso_idx', {}, 'abs_idx', {}, 'reason', {});

fprintf('Sampling logN values for each absorber...\n');
count = 0; skipped = 0;

for qso_idx = 1:size(logN_matrix, 1)
    logN_row = logN_matrix(qso_idx, :);
    p_row    = p_matrix(qso_idx, :);
    valid_absorbers = find(~isnan(logN_row) & logN_row > 0 & p_row > 0.99);

    for a = valid_absorbers
        if qso_idx > nqso_h5 || a > max_abs_per_qso
            skipped = skipped + 1;
            reasonsSkipped(skipped).qso_idx = qso_idx;
            reasonsSkipped(skipped).abs_idx = a;
            reasonsSkipped(skipped).reason = 'index out of bounds';
            continue;
        end

        try
            log_likelihoods = h5read(m.Properties.Source, ...
                '/LikelihoodCat/all_sample_log_likelihoods_c4L2', ...
                [qso_idx, 1, a], [1, nweights, 1]);
        catch
            skipped = skipped + 1;
            reasonsSkipped(skipped).qso_idx = qso_idx;
            reasonsSkipped(skipped).abs_idx = a;
            reasonsSkipped(skipped).reason = 'HDF5 read fail';
            continue;
        end

        log_likelihoods = double(log_likelihoods(:));

        if any(isnan(log_likelihoods)) || any(~isfinite(log_likelihoods))
            skipped = skipped + 1;
            reasonsSkipped(skipped).qso_idx = qso_idx;
            reasonsSkipped(skipped).abs_idx = a;
            reasonsSkipped(skipped).reason = 'invalid log-likelihoods';
            continue;
        end

        % Convert log-likelihoods to safe probability weights
        weights = exp(log_likelihoods - max(log_likelihoods));
        weights_sum = sum(weights);
        if weights_sum <= 0
            skipped = skipped + 1;
            reasonsSkipped(skipped).qso_idx = qso_idx;
            reasonsSkipped(skipped).abs_idx = a;
            reasonsSkipped(skipped).reason = 'zero weights after exp';
            continue;
        end
        weights = weights / weights_sum;

        try
            idx_samples = randsample(1:nweights, 1000, true, weights);
        catch
            skipped = skipped + 1;
            reasonsSkipped(skipped).qso_idx = qso_idx;
            reasonsSkipped(skipped).abs_idx = a;
            reasonsSkipped(skipped).reason = 'randsample fail';
            continue;
        end

        logN_samples = logN_grid(idx_samples);

        count = count + 1;
        sampledAbsorbers(count).qso_idx = qso_idx;
        sampledAbsorbers(count).abs_idx = a;
        sampledAbsorbers(count).logN_samples = logN_samples;

        if mod(count, 10000) == 0
            fprintf('Sampled %d absorbers (latest QSO %d)\n', count, qso_idx);
        end
    end

    if mod(qso_idx, 1000) == 0
        fprintf('Processed QSO %d/%d\n', qso_idx, size(logN_matrix, 1));
    end
end

%% Save results
save('sampled_logN_per_absorber.mat', 'sampledAbsorbers', 'reasonsSkipped', '-v7.3');
fprintf('\n Saved %d absorbers, %d skipped\n', count, skipped);