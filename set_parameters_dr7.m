% set_parameters: sets various parameters for the CIV detection 
% pipeline

% physical constants
civ_1548_wavelength = 1548.2049;		 % CIV transition wavelength  Å
civ_1550_wavelength =  1550.77845; 		 % CIV transition wavelength  Å
speed_of_light = 299792458;                   % speed of light                     m s⁻¹

% converts relative velocity in km s^-1 to redshift difference
kms_to_z = @(kms) (kms * 1000) / speed_of_light;

% utility functions for redshifting
emitted_wavelengths = ...
    @(observed_wavelengths, z) (observed_wavelengths / (1 + z));

observed_wavelengths = ...
    @(emitted_wavelengths,  z) ( emitted_wavelengths * (1 + z));

release = 'dr7';
% download Cooksey's dr7 spectra from this page: 
% http://www.guavanator.uhh.hawaii.edu/~kcooksey/SDSS/CIV/index.html 
% go to table: "SDSS spectra of the sightlines surveyed for C IV."
file_loader = @(mjd, plate, fiber_id) ...
  (read_spec_dr7(sprintf('data/dr7/spectro/1d_26/%04i/1d/spSpec-%05i-%04i-%03i.fit',...
  plate, mjd, plate, fiber_id)));


% file loading parameters
loading_min_lambda = 1315;          % range of rest wavelengths to load  Å
loading_max_lambda = 1550;                    
% The maximum allowed is set so that even if the peak is redshifted off the end, the
% quasar still has data in the range

% preprocessing parameters
%z_qso_cut      = 2.15;                   % filter out QSOs with z less than this threshold
z_qso_cut      = 1.7;                      % according to Cooksey z>1.7                      
min_num_pixels = 200;                         % minimum number of non-masked pixels

% normalization parameters
% range of rest wavelengths to use   Å
normalization_min_lambda = 1417; 
normalization_max_lambda = 1486; 
% null model parameters
min_lambda         =  1295;                   % range of rest wavelengths to       Å
max_lambda         = 1570;                    %   model
dlambda            = 0.9;                    % separation of wavelength grid      Å
k                  = 20;                      % rank of non-diagonal contribution
max_noise_variance = 1^2;                    % maximum pixel noise allowed during model training
masking_CIV_region = 0;                       % handel to mask CIV absorption to increase the training sample size
h                  = 2;                       % masking par to remove CIV region 
nAVG               = 0;                       % number of points added between two 
                                             % observed wavelengths to make the Voigt finer
SkyLine            = 0;                       % Handel for removing/not removing sky lines as indicated in C13                                             
                                        
% optimization parameters
minFunc_options =               ...           % optimization options for model fitting
    struct('MaxIter',     10000, ...
           'MaxFunEvals', 10000);

% C4 model parameters: parameter samples (for Quasi-Monte Carlo)
RejectionSampling        = 0;                      % Rejection Sampling handle 
num_C4_samples           = 50000;                  % number of parameter samples
alpha                    = 0.9;                    % weight of KDE component in mixture
uniform_min_log_nciv     = 14.0;                   % range of column density samples    [cm⁻²]
uniform_max_log_nciv     = 15.8;                   % from uniform distribution
fit_min_log_nciv         = 14.0;                   % range of column density samples    [cm⁻²]
fit_max_log_nciv         = 15.8;                   % from fit to log PDF
extrapolate_min_log_nciv = 14.0;  

min_sigma                = 20e5;                   % cm/s -> b/sqrt(2) -> min Doppler par from Cooksey
max_sigma                = 65e5;                   % cm/s -> b/sqrt(2) -> max Doppler par from Cooksey
% vCut                     = 3000;                    % maximum cut velocity for CIV system 
vCut                     = 5000;                    % maximum cut velocity for CIV system 
% model prior parameters
prior_z_qso_increase = kms_to_z(30000);       % use QSOs with z < (z_QSO + x) for prior

% instrumental broadening parameters
width = 3;                                    % width of Gaussian broadening (# pixels)
pixel_spacing = 1e-4; 

% wavelength spacing of pixels in dex

% CIV model parameters: absorber range and model
num_lines = 2;                                % number of members of CIV series to use

max_z_cut = kms_to_z(vCut);                   % max z_DLA = z_QSO - max_z_cut
% max_z_c4 = @(wavelengths, z_qso) ...         % determines maximum z_DLA to search
%     (max(wavelengths)/civ_1548_wavelength - 1) - max_z_cut;
max_z_c4 = @(z_qso, max_z_cut) ...         % determines maximum z_civ to search
     z_qso - max_z_cut*(1+z_qso);
min_z_cut = kms_to_z(vCut);                   % min z_DLA = z_Ly∞ + min_z_cut
% min_z_c4 = @(wavelengths, z_qso) ...         % determines minimum z_DLA to search
%     max(min(wavelengths) / civ_1548_wavelength - 1,                          ...
%         observed_wavelengths(min_lambda, z_qso) / civ_1548_wavelength - 1 + ...
%         min_z_cut);
min_z_c4 = @(wavelengths, z_qso) ...         % determines minimum z_civ to search
     min(wavelengths) / civ_1548_wavelength - 1;
train_ratio =0.95;
training_set_name =sprintf('lmin-%d-lmax-%d-Skyline-%d-norm-%d-%d-k-%d-dl-%d-civ-%d',... 
min_lambda, max_lambda, SkyLine, normalization_min_lambda,...
normalization_max_lambda, k, floor(dlambda*100), masking_CIV_region);
% base directory for all data
base_directory = 'data';
% utility functions for identifying various directories
distfiles_directory = @(release) ...
   sprintf('%s/%s/distfiles', base_directory, release);

spectra_directory   = @(release)...
   sprintf('%s/%s/spectra', base_directory, release);

processed_directory = @(release) ...
   sprintf('%s/%s/processed', base_directory, release);

c4_catalog_directory = @(name) ...
   sprintf('%s/C4_catalogs/%s/processed', base_directory, name);

   
% replace with @(varargin) (fprintf(varargin{:})) to show debug statements
% fprintf_debug = @(varargin) (fprintf(varargin{:}));
% fprintf_debug = @(varargin) ([]);
