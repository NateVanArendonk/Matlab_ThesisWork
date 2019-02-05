clearvars 
clc

% Load in Raster and polygon of AOI 
load('E:\Abbas\WCRP\Tacoma\Tier3\FloodRaster\RasterOUT\SWL\OBRW_Contemporary_SWL_FloodRaster.mat')
c = load('E:\Abbas\WCRP\Tacoma\Tier3\FloodRaster\SWL_PolygonMasks\Contemporary_polygon_SWL_RasterMask.mat');
%% Convet to lat and lon 
x = x(:);
y = y(:);

utmzone = repmat('10 T',size(x));
tic
[lat,lon] = utm2deg(x,y,utmzone);
toc

%% Reshape 
x = reshape(x,size(z));
y = reshape(y,size(z));

lat = reshape(lat,size(z));
lon = reshape(lon,size(z));
%% Subsample WL 

wl(~in) = NaN;
flood = wl - z;
flood(flood<0) = NaN;
%%
% Make blue colormap
cmap = colormap(jet(100));
cmap = cmap(90:-1:60,:,:);

cLimLow = 3;
cLimHigh = 4;
altitude = 10000;
alphaMatrix = ones(size(lon))*.55;

kmlFileName = 'test_geToolbox.kml';

output = ge_imagesc(lon(:),lat(:),flood,...
    'cLimLow',cLimLow,'cLimHigh',cLimHigh,...
    'altitude',altitude,...
    'altitudeMode','absolute',...
    'colorMap',cmap,...
    'alphaMatrix',alphaMatrix);

output2 = ge_colorbar(lon(end),lat(1),flood,...
    'numClasses',10,...
    'cLimLow',cLimLow,...
    'cLimHigh',cLimHigh,...
    'cBarFormatStr','%+0.1f',...
    'colorMap',cmap);

ge_output(kmlFileName,[output2 output],'name',kmlFileName);
