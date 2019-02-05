%% Build a Tide data Structure
%% Load in data
clearvars
stn_nm = 'seattle/';

dir_nm = '../COOPS_tides/';

fol_loc = strcat(dir_nm, stn_nm);
% Load predictions
pre = load(strcat(fol_loc,'seattle_tide_predictions'));

%Load verified tides
ver_t = load(strcat(fol_loc,'seattle_6minV'));

ver.time = ver_t.tides.time;
ver.wl = ver_t.tides.WL_VALUE;
ver.stn_id = ver_t.tides.STATION_ID;
ver.lat = ver_t.tides.latitude;
ver.lon = ver_t.tides.longitude;

clear ver_t



%%
% Change verified to ver
% Change raw to raw
% Change predicted to pr

tides = struct();

tides.lat = ver.lat;
tides.lon = ver.lon;
tides.datum = 'NAVD88';
tides.stn_id = ver.stn_id;

tides.ver_time = ver.time;
%tides.raw_time = raw.time;
tides.pre_time = pre.time;

tides.ver_wl = ver.wl;
%tides.raw_wl = raw.WL_VALUE;
tides.pre_wl = pre.wl;



%% Calculate Non-Tidal Residual

% First create a time vector and empty waterlevel vector
tic

t1 = pre.time(1);   % starting time
t2 = pre.time(end);   % ending time
t_vec = t1:(6/24/60):t2;   % time vector
temp_wl = NaN(length(t_vec), 1);   % water level empty vector

tic
% Put data onto vector with NaN'd gaps instead of missing data
for j = 1:10000%length(t_vec)
    I = find(t_vec(j)==ver.time);
    if ~isempty(I)
        temp_wl(j) = ver.wl(I(1));
    end
end
    
toc

%%
clf
plot(t_vec(1:10000), temp_wl(1:10000))
hold on
plot(pre.time(1:10000)-8/24, pre.wl(1:10000))