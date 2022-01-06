% learn_qso_model: fits GP to training catalog via maximum likelihood

rng('default');

% determine which spectra to use for training; allow string value for
% train_ind
if (ischar(train_ind))
  train_ind = eval(train_ind);
end

% select training vectors
all_wavelengths    =    all_wavelengths(train_ind, :);
all_wavelengths    =          all_wavelengths(mask_c4);
all_flux           =           all_flux(train_ind, :);
all_flux           =                all_flux(mask_c4);
all_noise_variance = all_noise_variance(train_ind, :);
all_noise_variance =       all_noise_variance(mask_c4);
all_pixel_mask     =     all_pixel_mask(train_ind, :);
all_pixel_mask     =           all_pixel_mask(mask_c4);
% z_qsos             =        all_zqso(train_ind);
z_qsos             = z_qsos(mask_c4);

num_quasars = numel(z_qsos);
fprintf('#QSOs %i\n', num_quasars);

rest_wavelengths = (min_lambda:dlambda:max_lambda);
num_rest_pixels  = numel(rest_wavelengths);

rest_fluxes          = nan(num_quasars, num_rest_pixels);
rest_noise_variances = nan(num_quasars, num_rest_pixels);
restSN               = nan(num_quasars,1);

% the preload_qsos should fliter out empty spectra;
% this line is to prevent there is any empty spectra
% in preloaded_qsos.mat for some reason
is_empty             = false(num_quasars, 1);

% interpolate quasars onto chosen rest wavelength grid
for i = 1:num_quasars
  z_qso = z_qsos(i);

  this_wavelengths    =    all_wavelengths{i}';
  this_flux           =           all_flux{i}';
  this_noise_variance = all_noise_variance{i}';
  this_pixel_mask     =     all_pixel_mask{i}';

  this_rest_wavelengths = emitted_wavelengths(this_wavelengths, z_qso);

  this_flux(this_pixel_mask)           = nan;
  this_noise_variance(this_pixel_mask) = nan;
  if (masking_CIV_region == 1)
    % masking CIV absorption 
    l1 = 1548.2040;
    l2 = 1550.7810;
    maskC4 = this_pixel_mask;
    nMasked=0;
    for Sys=1:17 
      if(all_z_civ(i,Sys)>0)
        maskC4 = this_wavelengths>(1+all_z_civ(i,Sys))*(l1-h) & this_wavelengths<(1+all_z_civ(i,Sys))*(l2+h);
        this_flux(maskC4)=nan;
        this_noise_variance(maskC4)=nan; 
        nMasked =nMasked+sum(maskC4);
        
      end
    end
      fprintf('processing quasar %i with lambda_size = %i %i, CIV-Masked:%d...\n', i, size(this_wavelengths), nMasked);
  else
    fprintf('processing quasar %i with lambda_size = %i %i \n', i, size(this_wavelengths));
  

  
  end
  
  
  if all(size(this_wavelengths) == [0 0])
    is_empty(i, 1) = 1;
    continue;
  end

  rest_fluxes(i, :) = ...
      interp1(this_rest_wavelengths, this_flux,           rest_wavelengths);

  % normalizing here
  ind = (this_rest_wavelengths >= normalization_min_lambda) & ...
          (this_rest_wavelengths <= normalization_max_lambda) & ...
          (~this_pixel_mask);
 
this_median = nanmedian(this_flux(ind));
  rest_fluxes(i, :) = rest_fluxes(i, :) / this_median;

  rest_noise_variances(i, :) = ...
      interp1(this_rest_wavelengths, this_noise_variance, rest_wavelengths);
  rest_noise_variances(i, :) = rest_noise_variances(i, :) / this_median .^ 2;

  restSN(i,1) = nanmedian(this_rest_wavelengths./sqrt(this_noise_variance));

  % Following C13 median signal to noise ratio in the region of CIV 
  % search should be larger than 4
  if restSN(i)>4
    filter_flags(i) = bitset(filter_flags(i), 5, true);
  end 
end
clear('all_wavelengths', 'all_flux', 'all_noise_variance', 'all_pixel_mask');

% filter out empty spectra
% note: if you've done this in preload_qsos then skip these lines
z_qsos               = z_qsos(~is_empty);
rest_fluxes          = rest_fluxes(~is_empty, :);
rest_noise_variances = rest_noise_variances(~is_empty, :);

% update num_quasars in consideration
num_quasars = numel(z_qsos);

fprintf('Get rid of %i empty spectra.\nnum_quasars = %i\n', sum(is_empty), num_quasars);

% mask noisy pixels
ind = (rest_noise_variances > max_noise_variance);

fprintf("Masking %g of pixels\n", nnz(ind) * 1 ./ numel(ind));

rest_fluxes(ind)          = nan;
rest_noise_variances(ind) = nan;




% Filter out spectra which have too many NaN pixels
ind = sum(isnan(rest_fluxes),2) < num_rest_pixels-min_num_pixels;

fprintf("Filtering %g quasars for NaN\n", num_quasars - nnz(ind));

rest_fluxes          = rest_fluxes(ind, :);
rest_noise_variances = rest_noise_variances(ind,:);

% Check for columns which contain only NaN on either end.
nancolfrac = sum(isnan(rest_fluxes), 1) / nnz(ind);

fprintf("Columns with nan > 0.9: ");

max(find(nancolfrac > 0.9))

% find empirical mean vector and center data
mu = nanmean(rest_fluxes);
centered_rest_fluxes = bsxfun(@minus, rest_fluxes, mu);

clear('rest_fluxes');

% small fix to the data fit into the pca:
% make the NaNs to the medians of a given row
% rememeber not to inject this into the actual
% joint likelihood maximisation
pca_centered_rest_flux = centered_rest_fluxes;

[num_quasars, ~] = size(pca_centered_rest_flux);

for i = 1:num_quasars
  this_pca_centered_rest_flux = pca_centered_rest_flux(i, :);

  % assign median value for each row to nan
  ind = isnan(this_pca_centered_rest_flux);

  pca_centered_rest_flux(i, ind) = nanmedian(this_pca_centered_rest_flux);
end

% get top-k PCA vectors to initialize M
[coefficients, ~, latent] = ...
  pca(pca_centered_rest_flux, ...
        'numcomponents', k, ...
        'rows',          'complete')

% initialize A to top-k PCA components of non-DLA-containing spectra
initial_M = bsxfun(@times, coefficients(:, 1:k), sqrt(latent(1:k))');
objective_function = @(x) objective(x, centered_rest_fluxes, rest_noise_variances);

% maximize likelihood via L-BFGS
[x, log_likelihood, ~, minFunc_output] = ...
    minFunc(objective_function, initial_M, minFunc_options);

ind = (1:(num_rest_pixels * k));
M = reshape(x(ind), [num_rest_pixels, k]);

variables_to_save = {'release', 'train_ind', 'max_noise_variance', ...
                     'minFunc_options', 'rest_wavelengths', 'mu', ...
                     'initial_M', 'M',  'log_likelihood', ...
                     'minFunc_output',   'restSN', 'coefficients',...
                     'latent'};

save(sprintf('%s/learned_model-%s', processed_directory(release), ...
                                      training_set_name),...
            variables_to_save{:}, '-v7.3');
