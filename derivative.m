% Total Horizontal Derivative (THD) 
%
% This script computes spatial gravity derivatives in the frequency domain 
% from regional and residual Bouguer anomaly data, including:
%   - Gzx: Horizontal derivative in x-direction
%   - Gzy: Horizontal derivative in y-direction
%   - Gzz: Vertical derivative
%   - THD: Total Horizontal Derivative (√(Gzx² + Gzy²))
%
% This script performs following steps:
% 1. Reads the regional or residual gravity field (GeoTIFF)
% 2. Performs FFT and applies a Gaussian low-pass filter to reduce noise
% 3. Calculates derivatives in the frequency domain using spectral operators
% 4. Inverse transforms to recover spatial domain fields
% 5. Crops edge artifacts caused by FFT wrapping
% 6. Outputs and visualizes the THD and derivative fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
%% import gravity data (Regional_CBA and Residual_CBA)
[G, R] = readgeoraster('data\Regional_CBA.tif'); 

%% perform Fourier transform
G_fft = inpaint_nans(G, 4);

G_fft = fftshift(fft2(G_fft));

%% low pass filter to eliminate noise 
[X, Y] = meshgrid(-size(G_fft,2)/2: size(G_fft,2)/2-1, -size(G_fft,1)/2: size(G_fft,1)/2-1);
sigma = 40;
low_pass = exp(-(X.^2 + Y.^2) / (2 * sigma^2));

G_fft = G_fft .* low_pass;
%% Calculate wavenumber
[kx, ky] = meshgrid((0:size(G_fft,2)-1) - floor(size(G_fft,2)/2), ...
                    (0:size(G_fft,1)-1) - floor(size(G_fft,1)/2)); %frequencies indices
% note : GGM plus has resolution of 7.2 arc second = 220 meters at 8.3 deg
% calculate project lenght
kx = kx.*(2*pi/(R.RasterSize(2)*220));
ky = ky.*(2*pi/(R.RasterSize(1)*220));
%kx = kx .* (2*pi/ range(xq));
%ky = ky .* (2*pi/ range(yq));
k = sqrt(kx.^2 + ky.^2);

%% xyz derivative
Gzx = 1i.*kx.*G_fft;
Gzy = 1i.*ky.*G_fft;
Gzz = k .* G_fft;

%% inverse fourier to spatial domain 
Gzx = real(ifft2(ifftshift(Gzx)));
Gzy = real(ifft2(ifftshift(Gzy)));
Gzz = real(ifft2(ifftshift(Gzz)));

Gzx(isnan(G)) = NaN;
Gzy(isnan(G)) = NaN;
Gzz(isnan(G)) = NaN;

%% crop data (crop out the artifact)
crop = 0.05; %%crop percent
Gzx = Gzx(ceil(crop*size(Gzx,1)) : end - floor(crop*size(Gzx,1)), ...
            ceil(crop*size(Gzx,2)) : end - floor(crop*size(Gzx,2)));

Gzy = Gzy(ceil(crop*size(Gzy,1)) : end - floor(crop*size(Gzy,1)), ...
            ceil(crop*size(Gzy,2)) : end - floor(crop*size(Gzy,2)));

Gzz = Gzz(ceil(crop*size(Gzz,1)) : end - floor(crop*size(Gzz,1)), ...
            ceil(crop*size(Gzz,2)) : end - floor(crop*size(Gzz,2)));

%% create new reference 
latlim = [R.LatitudeLimits(1)+(crop*range(R.LatitudeLimits)) , R.LatitudeLimits(2)-(crop*range(R.LatitudeLimits))];
lonlim = [R.LongitudeLimits(1)+(crop*range(R.LongitudeLimits)) ,  R.LongitudeLimits(2)-(crop*range(R.LongitudeLimits))];
rasterSize = [size(Gzx,1), size(Gzx,2)];
R_new = georefcells(latlim,lonlim,rasterSize,'ColumnsStartFrom','north');

%% plot
figure;
subplot 131 ;
geoshow(Gzx,R_new,'DisplayType','surface');
colorbar;
subplot 132 ;
geoshow(Gzy,R_new,'DisplayType','surface');
colorbar;
subplot 133;
geoshow(Gzz,R_new,'DisplayType','surface');
colorbar;
colormap turbo

figure;
THD = sqrt(Gzx.^2 + Gzy.^2);
geoshow(THD, R_new, 'DisplayType','surface');
colormap turbo
colorbar

%% Export GeoTiff
geotiffwrite('data/Regional_THD.tif', THD, R_new); %Change the output file name as necessary.