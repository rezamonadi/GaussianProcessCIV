close all 
fig = figure()
sigma =10e5;
map_z = 3.75;
map_N= 13.5;
padded_wavelengths = ...
    [logspace(log10(min(this_unmasked_wavelengths)) - width * pixel_spacing, ...
    log10(min(this_unmasked_wavelengths)) - pixel_spacing,...
    width)';...
    this_unmasked_wavelengths;...
    logspace(log10(max(this_unmasked_wavelengths)) + pixel_spacing,...
    log10(max(this_unmasked_wavelengths)) + width * pixel_spacing,...
    width)'...
    ];
wv = logspace(log10(min(this_unmasked_wavelengths)), log10(max(this_unmasked_wavelengths)), 10000);

padded_wv = ...
    [logspace(log10(min(this_unmasked_wavelengths)) - width * pixel_spacing, ...
    log10(min(this_unmasked_wavelengths)) - pixel_spacing,...
    width)';...
    wv';...
    logspace(log10(max(this_unmasked_wavelengths)) + pixel_spacing,...
    log10(max(this_unmasked_wavelengths)) + width * pixel_spacing,...
    width)'...
    ];

padded_sigma_pixels = ...
        [this_sigma_pixel(1)*ones(width,1);...
        this_sigma_pixel;...
        this_sigma_pixel(end)*ones(width,1)];

a_fine = voigt_noB(padded_wv, map_z, ...
10^map_N,num_lines, sigma);
a_org = voigt_noB(padded_wavelengths, map_z, ...
10^map_N,num_lines, sigma);

plot(wv, a_fine, 'LineWidth', 1.5, 'Marker', 'o')
hold on
plot(this_unmasked_wavelengths, a_org, 'Marker', 'x', 'LineWidth',1.5)
hold on
for nAVG = [2,  20]
    padded_wavelengths_fine = ...
    [logspace(log10(min(this_unmasked_wavelengths)) - width * pixel_spacing/(nAVG+1), ...
    log10(min(this_unmasked_wavelengths)) - pixel_spacing/(nAVG+1),...
    width)';...
    finer(this_unmasked_wavelengths, nAVG)';...
    logspace(log10(max(this_unmasked_wavelengths)) + pixel_spacing/(nAVG+1),...
    log10(max(this_unmasked_wavelengths)) + width * pixel_spacing/(nAVG+1),...
    width)'...
    ];

    padded_sigma_pixels_fine = ...
        [this_sigma_pixel(1)*ones(width,1);...
        finer(this_sigma_pixel, nAVG)';...
        this_sigma_pixel(end)*ones(width,1)];

    
    a_fine = voigt_noB(padded_wavelengths_fine, map_z, ...
                10^map_N,num_lines, sigma);

    a_avg = Averager(a_fine, nAVG, lenW_unmasked);
    z_civ= this_unmasked_wavelengths/civ_1548_wavelength -1; 
    % indPlot = abs(z_civ - map_z)<0.00001;
    p= plot(this_unmasked_wavelengths, a_avg);
    % nnz(indPlot)
    p.LineWidth=1.5;
    p.Marker = 'x';
    hold on 
end
legend('fine', 'nAVG=0','nAVG=2', 'nAVG=20', 'location','southwest');
xlim([7345, 7375])

xlabel('$\lambda(\AA)$', 'interpreter', 'latex')
ylabel('$a_{1548\AA}\times a_{1550\AA}$', 'interpreter', 'latex')
exportgraphics(fig, 'VoigtTest_noB.png', 'Resolution', 800);