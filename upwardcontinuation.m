% Regional Anomaly Extraction via Upward Continuation
%
% This script performs upward continuation on a complete Bouguer anomaly
% using Fast Fourier Transform (FFT). The goal is to estimate the
% regional gravity field by attenuating shorter-wavelength.
%
% This script perform the following step
% 1. The complete Bouguer anomaly is imported from a GeoTIFF
% 2. 2D Fourier transform is applied to shift to frequency domain
% 3. A low-pass exponential filter (based on Blakely, 1995) is used
% 4. The result is transformed back to spatial domain to yield regional field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear 
clc

%% import gravity data 
[CBA, R] = readgeoraster('data/CompleteBouguerAnomaly.tif'); 

%% perform Fourier transform
G_fft = inpaint_nans(CBA, 4);
G_fft = fftshift(fft2(G_fft));

%% Define grid
xq = linspace(R.LongitudeLimits(1), R.LongitudeLimits(2), R.RasterSize(2));
yq = linspace(R.LatitudeLimits(1), R.LatitudeLimits(2), R.RasterSize(1));
[Xq, Yq] = meshgrid(xq, yq);

%% Calculate wavenumber
[kx, ky] = meshgrid((0:size(Xq,2)-1) - floor(size(Xq,2)/2), ...
                    (0:size(Yq,1)-1) - floor(size(Yq,1)/2)); %frequencies indices
% note : GGM plus has resolution of 7.2 arc second = 220 meters at 8.3 deg
% calculate project lenght
kx = kx.*(2*pi/(R.RasterSize(2)*220));
ky = ky.*(2*pi/(R.RasterSize(1)*220));
k = sqrt(kx.^2 + ky.^2);

%% Upward Continuation (Blakely, 1995)
h = 5000; 
H = exp(-h*k); %Operator
G_reg = G_fft .* H;
G_reg = ifftshift(G_reg); 

G_reg = real(ifft2(G_reg));
G_reg(isnan(CBA)) = NaN; % clip

%% Plot result
figure;
subplot 121 ;
geoshow(CBA, R, 'DisplayType','surface');
title('Bouguer Anomaly');
colorbar;

subplot 122 ;
geoshow(G_reg, R, 'DisplayType','surface');
title('Reigonal Anomaly');
colorbar;
colormap turbo;

%% Export 
geotiffwrite('data\Regional_CBA.tif', G_reg, R)