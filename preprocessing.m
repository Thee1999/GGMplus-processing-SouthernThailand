% Gravity Preprocessing Script for Bouguer Anomaly Calculation
%
% This script performs the following steps
% 1. Import Free Air anomaly and DEM 
% 2. Compute the Simple Bouguer anomaly 
% 3. Applies terrain correction derived from terraincorrection.m
% 4. Plots all product and export the Complete Bouguer Anomaly as GeoTiff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear

%% import Free Air anomaly
[FAA, R1]= readgeoraster('data/FreeAirAnomaly.tif'); % Free Air anomaly derived from GGMplus data
[DEM, R2]= readgeoraster('data/demGGM.tif'); % Import DEM with identical grid with GGMplus model
G = 6.6742e-11;%Gravity constant 

%% Bouguer Correction
% Calculate bimple Bouguer anomaly assuming a density of Bouguer slab of
% 2.67 g/cm3
% The Bouguer correction was calculated as Δg_B = 2πGρh simplified to Δg_B ≈ 0.1119 × h (mGal)

SBA = FAA;
SBA = SBA - 0.1119.*DEM; % Simple Bouguer anomaly (mGal) 

%% Terrain Correction
[TC, R3] = readgeoraster('data/TerrainCorrection.tif');
CBA = SBA;
CBA = CBA + TC;

%% Plot
figure;

subplot 131;
geoshow(FAA, R1, 'DisplayType','surface');
title('Free Air Anomaly');
colorbar ; 

subplot 132;
geoshow(SBA, R2, 'DisplayType','surface');
title('Simple Bouguer Anomaly');
colorbar ; 

subplot 133;
geoshow(CBA, R3, 'DisplayType','surface');
title('Complete Bouguer Anomaly')
colorbar ;
colormap 'turbo';
%% Export Geotiff
geotiffwrite('data/CompleteBouguerAnomalytest.tif', CBA, R3);