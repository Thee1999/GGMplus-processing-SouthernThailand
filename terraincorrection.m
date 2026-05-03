% Terrain Correction for Gravity Data Using DEM and Station Coordinates
%
% This script computes the terrain correction (TC) for each gravity station 
% using a high-resolution Digital Elevation Model (DEM) and a multi-zone 
% correction method based on wedge, prism, and sector approximations.
%
% Methods:
% - Inner Zone (180 m radius): Wedge approximation (Campbell, 1980)
% - Zones 1–2 (180 m to 2880 m): Prism approximation (Nagy, 1966)
% - Zones 3–7 (2880 m to 92160 m): Concentric square segments (Kane, 1962)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
clc 
format longg 
%% import DEM  // resolution 90 m 
DEM = load('DEM.mat');
X = gpuArray(DEM.terrain{1});
Y = gpuArray(DEM.terrain{2});
Z = gpuArray(DEM.terrain{3});

%% import station 
station = load("station.mat");
st = gpuArray(station.station_points(:,:));

%% Terrain Correction
terrain_effect = gpuArray(zeros(size(st,1),1));

for i = 2337501:size(st,1)
   terrain_effect(i,:) = TC(st(i,:),X,Y,Z) ; 
end

%% Define function terrain correction 
function TerrainCorr = TC(st,X,Y,Z)  
    %% Define Constant Parameter 
    density = 2670 ; %kg/m3
    G = 6.6742e-11 ; % Gravitational constant in m3/kg.s2
    si_to_mgal = 1e5 ; %conversion factor to mgal 

    if ~isequal(size(st), [1, 2])
       error('Input must contains x and y');
    end
%% locate station point
    X_AOI = X(1840:8560, 1811:7880);
    Y_AOI = Y(1840:8560, 1811:7880);
    distances = sqrt((X_AOI - st(1,1)).^2 + (Y_AOI - st(1,2)).^2);
    [row, col] = find(distances == min(distances,[] ,'all'));
    row = row +1839;
    col = col +1810;
%% Inner Zone within 2 cell, 180 m  (Campbell, 1980)
    Inner = Z(row-2:row+2, col-2:col+2);
    
    wedge(1) = mean(Inner(1:2, end-1:end),'all') - Z(row,col);
    wedge(2) = mean(Inner(end-1:end, end-1:end),'all') - Z(row,col);
    wedge(3) = mean(Inner(end-1:end, 1:2),'all') - Z(row,col);
    wedge(4) = mean(Inner(1:2, 1:2),'all') - Z(row,col); 
    
    theta = atand(wedge./180);
    gi = 2*pi*G*density*180.*(1-cosd(theta)).*si_to_mgal;
    TerrainCorr = sum(gi, 'all');
    clear Inner wedge theta gi
%% Zone 1 inner ring 180 m outer ring 1440 m, 16 row adjacent to the cell (Nagy, 1966)
    dh= Z(row-16:row+16, col-16:col+16) - Z(row,col);
    dX = X(row-16:row+16, col-16:col+16) - X(row, col);
    dY = Y(row-16:row+16, col-16:col+16) - Y(row, col);
    dh(ceil(end/2)-2:ceil(end/2)+2, ceil(end/2)-2:ceil(end/2)+2) = NaN;
    dh(17:17,:) = []; dh(:,17:17) = [];
    dX(17:17,:) = []; dX(:,17:17) = [];
    dY(17:17,:) = []; dY(:,17:17) = [];
    x1 = dX(1:2:end,2:2:end); x2 = dX(2:2:end, 1:2:end);
    y1 = dY(1:2:end,2:2:end); y2 = dY(2:2:end, 1:2:end);
    dh = padarray(dh, [1 1], 0,'both') ;
    dhmean = zeros(16,16);
    %calculate mean
    for i = 2:2:size(dh,1)-1
        for j = 2:2:size(dh,2)-1
            dhmean(floor(i/2),floor(j/2)) = sum(dh(i:i+1,j:j+1), "all", 'includenan')./4;
        end
    end
    %% Nagy TC %correct the variable 
    fac11=sqrt(x1.^2+y1.^2);
    fac11h=sqrt(x1.^2+y1.^2+dhmean.^2);
    fac12=sqrt(x1.^2+y2.^2);
    fac12h=sqrt(x1.^2+y2.^2+dhmean.^2);
    fac21=sqrt(x2.^2+y1.^2);
    fac21h=sqrt(x2.^2+y1.^2+dhmean.^2);
    fac22=sqrt(x2.^2+y2.^2);
    fac22h=sqrt(x2.^2+y2.^2+dhmean.^2);
    fac2h=sqrt(y2.^2+dhmean.^2);
    fac1h=sqrt(y1.^2+dhmean.^2);
    y2h=y2.^2+dhmean.^2;
    y1h=y1.^2+dhmean.^2;
    terrc=x2.*(log( (y2+fac22)./(y2+fac22h))-...
        log( (y1+fac21)./(y1+fac21h))) -...
        x1.*( log( (y2+fac12)./(y2+fac12h))-log( (y1+fac11)./(y1+fac11h))) +...
        y2.*(log( (x2+fac22)./(x2+fac22h))-log( (x1+fac12)./(x1+fac12h))) -...
        y1.*(log( (x2+fac21)./(x2+fac21h))-log( (x1+fac11)./(x1+fac11h))) +...
        dhmean.*(asin( complex((y2h +y2.*fac22h)./( (y2+fac22h).*fac2h)))-...
        asin( complex((y2h +y2.*fac12h)./( (y2+fac12h).*fac2h)))-...
        asin( complex((y1h +y1.*fac21h)./( (y1+fac21h).*fac1h)))+...
        asin( complex((y1h +y1.*fac11h)./( (y1+fac11h).*fac1h))));
    terrc = 2.*pi.*G.*density.*terrc.*si_to_mgal;
    TerrainCorr = TerrainCorr + sum(terrc ,'all', 'omitnan');

    clear tercc dh dhmean dX dY x1 x2 y1 y2 fac11 fac11h fac12 fac12h fac21 fac21h fac22 fac22h fac2h fac1h y1h y2h 
    %% Zone 2 inner ring 1440, outer ring 2880 (Nagy, 1966)
    dh = Z(row-32:row+32, col-32:col+32) - Z(row,col);
    dX = X(row-32:row+32, col-32:col+32) - X(row, col);
    dY = Y(row-32:row+32, col-32:col+32) - Y(row, col);
    dh(ceil(end/2)-16:ceil(end/2)+16, ceil(end/2)-16:ceil(end/2)+16) = NaN;
    dh(33:33,:) = []; dh(:,33:33) = [];
    dX(33:33,:) = []; dX(:,33:33) = [];
    dY(33:33,:) = []; dY(:,33:33) = [];
    x1 = dX(1:4:end,4:4:end); x2 = dX(4:4:end, 4:4:end);
    y1 = dY(1:4:end,4:4:end); y2 = dY(4:4:end, 4:4:end);
    dh = padarray(dh, [3 3], 0,'both') ;
    dhmean = zeros(16,16);
    for i = 4:4:size(dh,1)-3
        for j = 4:4:size(dh,2)-3
            dhmean(floor(i/4),floor(j/4)) = sum(dh(i:i+3,j:j+3),'all','includenan')./16;
        end
    end
    %% Nagy TC %correct the variable 
    fac11=sqrt(x1.^2+y1.^2);
    fac11h=sqrt(x1.^2+y1.^2+dhmean.^2);
    fac12=sqrt(x1.^2+y2.^2);
    fac12h=sqrt(x1.^2+y2.^2+dhmean.^2);
    fac21=sqrt(x2.^2+y1.^2);
    fac21h=sqrt(x2.^2+y1.^2+dhmean.^2);
    fac22=sqrt(x2.^2+y2.^2);
    fac22h=sqrt(x2.^2+y2.^2+dhmean.^2);
    fac2h=sqrt(y2.^2+dhmean.^2);
    fac1h=sqrt(y1.^2+dhmean.^2);
    y2h=y2.^2+dhmean.^2;
    y1h=y1.^2+dhmean.^2;
    terrc=x2.*(log( (y2+fac22)./(y2+fac22h))-...
        log( (y1+fac21)./(y1+fac21h))) -...
        x1.*( log( (y2+fac12)./(y2+fac12h))-log( (y1+fac11)./(y1+fac11h))) +...
        y2.*(log( (x2+fac22)./(x2+fac22h))-log( (x1+fac12)./(x1+fac12h))) -...
        y1.*(log( (x2+fac21)./(x2+fac21h))-log( (x1+fac11)./(x1+fac11h))) +...
        dhmean.*(asin( complex((y2h +y2.*fac22h)./( (y2+fac22h).*fac2h)))-...
        asin( complex((y2h +y2.*fac12h)./( (y2+fac12h).*fac2h)))-...
        asin( complex((y1h +y1.*fac21h)./( (y1+fac21h).*fac1h)))+...
        asin( complex((y1h +y1.*fac11h)./( (y1+fac11h).*fac1h))));
    terrc = 2.*pi.*G.*density.*terrc.*si_to_mgal;
    TerrainCorr = TerrainCorr + sum(terrc, 'all', 'omitnan');
    clear tercc dh dhmean dX dY x1 x2 y1 y2 fac11 fac11h fac12 fac12h fac21 fac21h fac22 fac22h fac2h fac1h y1h y2h 

    %% Zone 3 square segment ring inner rad = 2880, outer rad = 5760 (Kane, 1962)
    dh = Z(row-64:row+64, col-64:col+64) - Z(row,col);
    dh(ceil(end/2)-32:ceil(end/2)+32, ceil(end/2)-32:ceil(end/2)+32) = NaN;
    dh(65:65,:) = []; dh(:,65:65) = [];
    dh = padarray(dh, [31 31], 0,'both') ;
    dhmean = zeros(4,4);
    for i = 32:32:size(dh,1)-31
        for j = 32:32:size(dh,2)-31
            dhmean(floor(i/32),floor(j/32)) = mean(dh(i:i+31,j:j+31),'all');
        end
    end
    terrc = G*density*(2880^2).*((2880+sqrt((2880^2)+dhmean.^2)-sqrt((5760^2)+dhmean.^2))./(5760^2 - 2880^2)).*si_to_mgal;
    TerrainCorr = TerrainCorr + sum(terrc, 'all', 'omitnan');

    %% Zone 4 square segment ring inner rad = 5760, outer rad = 11520 (Kane, 1962)
    R1 = 5760;
    R2 = 11520;
    dh = Z(row-128:row+128, col-128:col+128) - Z(row,col);
    dh(ceil(end/2)-64:ceil(end/2)+64, ceil(end/2)-64:ceil(end/2)+64) = NaN;
    dh(129:129,:) = []; dh(:,129:129) = [];
    dh = padarray(dh, [63 63], 0,'both') ;
    dhmean = zeros(4,4);
    for i = 64:64:size(dh,1)-63
        for j = 64:64:size(dh,2)-63
            dhmean(floor(i/64),floor(j/64)) = mean(dh(i:i+63,j:j+63),'all');
        end
    end
    terrc = G*density*((R2-R1)^2).*(((R2-R1)+sqrt((R1^2)+dhmean.^2)-sqrt((R2^2)+dhmean.^2))./(R2^2 - R1^2)).*si_to_mgal;
    TerrainCorr = TerrainCorr + sum(terrc, 'all', 'omitnan');
    %% Zone 5 square segment ring inner rad = 11520, outer rad = 23040 (Kane, 1962)
    R1 = 11520;
    R2 = 23040;
    dh = Z(row-256:row+256, col-256:col+256) - Z(row,col);
    dh(ceil(end/2)-128:ceil(end/2)+128, ceil(end/2)-128:ceil(end/2)+128) = NaN;
    dh(257:257,:) = []; dh(:,257:257) = [];
    dh = padarray(dh, [127 127], 0,'both') ;
    dhmean = zeros(4,4);
    for i = 128:128:size(dh,1)-127
        for j = 128:128:size(dh,2)-127
            dhmean(floor(i/128),floor(j/128)) = mean(dh(i:i+127,j:j+127),'all');
        end
    end
    terrc = G*density*((R2-R1)^2).*(((R2-R1)+sqrt((R1^2)+dhmean.^2)-sqrt((R2^2)+dhmean.^2))./(R2^2 - R1^2)).*si_to_mgal;
    TerrainCorr = TerrainCorr + sum(terrc, 'all', 'omitnan');

    %% Zone 6 square segment ring inner rad = 23040, outer rad = 46080 (Kane, 1962)
    R1 = 23040;
    R2 = 46080;
    dh = Z(row-512:row+512, col-512:col+512) - Z(row,col);
    dh(ceil(end/2)-256:ceil(end/2)+256, ceil(end/2)-256:ceil(end/2)+256) = NaN;
    dh(257:257,:) = []; dh(:,257:257) = [];
    dh = padarray(dh, [255 255], 0,'both') ;
    dhmean = zeros(4,4);
    for i = 256:256:size(dh,1)-255
        for j = 256:256:size(dh,2)-255
            dhmean(floor(i/256),floor(j/256)) = mean(dh(i:i+255,j:j+255),'all');
        end
    end
    terrc = G*density*((R2-R1)^2).*(((R2-R1)+sqrt((R1^2)+dhmean.^2)-sqrt((R2^2)+dhmean.^2))./(R2^2 - R1^2)).*si_to_mgal;
    TerrainCorr = TerrainCorr + sum(terrc, 'all', 'omitnan');

    %% Zone 7 square segment ring inner rad = 46080, outer rad = 92160 (Kane, 1962)
    R1 = 46080;
    R2 = 92160;
    dh = Z(row-1024:row+1024, col-1024:col+1024) - Z(row,col);
    dh(ceil(end/2)-512:ceil(end/2)+512, ceil(end/2)-512:ceil(end/2)+512) = NaN;
    dh(513:513,:) = []; dh(:,513:513) = [];
    dh = padarray(dh, [511 511], 0,'both') ;
    dhmean = zeros(4,4);
    for i = 512:512:size(dh,1)-511
        for j = 512:512:size(dh,2)-511
            dhmean(floor(i/512),floor(j/512)) = mean(dh(i:i+511,j:j+511),'all');
        end
    end
    terrc = G*density*((R2-R1)^2).*(((R2-R1)+sqrt((R1^2)+dhmean.^2)-sqrt((R2^2)+dhmean.^2))./(R2^2 - R1^2)).*si_to_mgal;
    TerrainCorr = TerrainCorr + sum(terrc, 'all', 'omitnan');
end

