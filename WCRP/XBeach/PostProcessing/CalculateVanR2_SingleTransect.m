% This code will calculate a R2 value at each transect given a TWL and wave
% hindcast
addpath C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\RunUp\RunupFormulas

clearvars
% ------------------- Load in Bathy/DEM
D = load('E:\Abbas\Modeling Resources\PS_DEM\Ruston_Way\RustonWayCONED_DEM.mat');

% ------------------- Load in transects
kml_fol = 'C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\KML\RunupTransects\OwenBeach\';
kml_nm = 'OB1.kml';
temp = kml2struct([kml_fol kml_nm]); % load in KML of Transect
T.lat = temp.Lat;
T.lon = temp.Lon;
[T.x_utm,T.y_utm] = deg2utm(T.lat,T.lon); % Convert to utm
T.name = kml_nm;
clear kml_fol kmls kml_nm


% ------------------- Load in wave hindcast
% W = load('C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\Hindcast\Ruston_LUT\OwenRuston_wave_hindcast_SLR1.3FT.mat');
W = load('C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\Hindcast\Ruston_LUT\OwenRuston_wave_hindcast.mat');
% Load in tides to get NTR
tide_fol = 'E:\Abbas\Model_Met_Forcings\COOPS_tides\';
tide_stns = 'tacoma\Tacoma_Reconstruct_NAVD88.mat';
tides = load([tide_fol tide_stns]);
T.twl = tides.twl_patched;
T.time = tides.time;
T.ntr = tides.ntr;

% Interpolate twl onto time vec
inds = find(unique(T.time));
T.time = T.time(inds);
T.twl = T.twl(inds);
T.ntr = T.ntr(inds);
inds = find(diff(T.time) == 0);
T.time(inds) = [];
T.twl(inds) = [];
T.ntr(inds) = [];
T.twl_i = double(interp1(T.time,T.twl,W.time));
T.ntr_i = double(interp1(T.time,T.ntr,W.time));

%% Calculate R2

V2 = zeros(1,length(W.time));
tic
% Make finely spaced transect from ends of transect
step = 1;
[ line_x, line_y ] = createTransect( T.x_utm(1), T.y_utm(1), T.x_utm(2), T.y_utm(2), step );
% Subset DEM to get elevations along transect
t_inds = D.x >= min(line_x) & D.x <= max(line_x) & ...
    D.y >= min(line_y) & D.y <= max(line_y);
t_x = D.x(t_inds);
t_y = D.y(t_inds);
t_z = D.z(t_inds);
% Get Elevations at each point
line_z = zeros(size(line_x)); % DEM Line
for jj = 1:length(line_x)
    dist = sqrt((line_x(jj)-t_x).^2 + (line_y(jj)-t_y).^2);
    [~, I] = min(dist);
    line_z(jj) = t_z(I);
end
% create along transect variable
s = sqrt(line_x.^2+line_y.^2);
s = s - min(s);
if s(1) ~= 0
    s = fliplr(s);
end
% Refine s and z to grab slope over small distances
s_fine = s(1):.0001:s(end);
z_fine = interp1(s,line_z,s_fine);

% Find Wave point that is closest to the transect
dist = sqrt((line_x(1)-W.x).^2 + (line_y(1)-W.y).^2);
[~, I] = min(dist);
% Get waves and period
hs = W.hs_ts(I,:);
tp = W.tp_ts(I,:);

% Calculate R2
for ii = 1:length(hs)
    V2(ii) = VanderMeer2002_R2(hs(ii),tp(ii),.15,-5,W.twl(ii)',0.80,s_fine,z_fine);
    if rem(ii,100) == 0
        fprintf('Completed %d out of %d\n',ii,length(hs))
    end
    
end

toc

%% Save
ntr = T.ntr_i;
twl = T.twl_i;
spd = W.speed;
wnddir = W.wnddir;

save('OB1_R2_rough80','V2','hs','tp','ntr','twl','spd','wnddir')%,'slr_val','slr_units')