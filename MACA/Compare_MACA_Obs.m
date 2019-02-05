% This script will load in MACA and Multiple obs stations and 
% 1. compare distributions for nearest MACA point to obs. to see if historic
% distributions really captures historic conditions.  
% 2. It will also plot The MACA points and the Obs points
% 3. Compare extremes with Maca and Obs


%% Load in Data 

% MACA First -------------------------------------------
% ------------------------ Load in MACA -----------------------------------
addpath C:\Functions_Matlab\time
addpath C:\Functions_Matlab
clearvars
% ---------------------- Historic Data ------------------------------------
Mh_u = ncread('C:\Users\ahooshmand\Desktop\Data_Forcings\MACA\maca_u_GFDL_ESM2M_historical_1950_2005.nc','eastward_wind');
Mh_v = ncread('C:\Users\ahooshmand\Desktop\Data_Forcings\MACA\maca_v_GFDL_ESM2M_historical_1950_2005.nc','northward_wind');
Mh_time = ncread('C:\Users\ahooshmand\Desktop\Data_Forcings\MACA\maca_u_GFDL_ESM2M_historical_1950_2005.nc','time');
Mh_time = double(Mh_time); % convert to double for following step
Mh_time = datenum(1900,1,1+Mh_time); % Convert MACA time to Matlab Datenum
Mh_lat = ncread('C:\Users\ahooshmand\Desktop\Data_Forcings\MACA\maca_u_GFDL_ESM2M_historical_1950_2005.nc','lat');
Mh_lon = ncread('C:\Users\ahooshmand\Desktop\Data_Forcings\MACA\maca_u_GFDL_ESM2M_historical_1950_2005.nc','lon');
Mh_lon = -1*(360 - Mh_lon);
Mh_spd = hypot(Mh_u,Mh_v);
wnddir_temp = (180/pi)*atan2(Mh_u,Mh_v);
% Rotate winds to compass directions
wnddir_temp = 90 - wnddir_temp;
wnddir_temp(wnddir_temp<0)=wnddir_temp(wnddir_temp<0)+360;
% switch to conventional coming from dir
wnddir_temp = wnddir_temp+180;
Mh_dir = wrap2360(wnddir_temp);
clear wnddir_temp

% ---------------------- Future Data --------------------------------------
scenario = 'rcp45';
switch scenario
    case 'rcp45'
        Mf_u = ncread('C:\Users\ahooshmand\Desktop\Data_Forcings\MACA\maca_u_GFDL_ESM2M_rcp45_2006_2100.nc','eastward_wind');
        Mf_v = ncread('C:\Users\ahooshmand\Desktop\Data_Forcings\MACA\maca_v_GFDL_ESM2M_rcp45_2006_2100.nc','northward_wind');
        Mf_time = ncread('C:\Users\ahooshmand\Desktop\Data_Forcings\MACA\maca_u_GFDL_ESM2M_rcp45_2006_2100.nc','time');
        Mf_time = double(Mf_time); % convert to double for following step
        Mf_time = datenum(1900,1,1+Mf_time); % Convert MACA time to Matlab Datenum
        Mf_spd = hypot(Mf_u,Mf_v);
        wnddir_temp = (180/pi)*atan2(Mf_u,Mf_v);
        % Rotate winds to compass directions
        wnddir_temp = 90 - wnddir_temp;
        wnddir_temp(wnddir_temp<0)=wnddir_temp(wnddir_temp<0)+360;
        % switch to conventional coming from dir
        wnddir_temp = wnddir_temp+180;
        Mf_dir = wrap2360(wnddir_temp);
        clear wnddir_temp
    case 'rcp85'
        Mf_u = ncread('C:\Users\ahooshmand\Desktop\Data_Forcings\MACA\maca_u_GFDL_ESM2M_rcp45_2006_2100.nc','eastward_wind');
        Mf_v = ncread('C:\Users\ahooshmand\Desktop\Data_Forcings\MACA\maca_v_GFDL_ESM2M_rcp45_2006_2100.nc','northward_wind');
        Mf_time = ncread('C:\Users\ahooshmand\Desktop\Data_Forcings\MACA\maca_u_GFDL_ESM2M_rcp45_2006_2100.nc','time');
        Mf_time = double(Mf_time); % convert to double for following step
        Mf_time = datenum(1900,1,1+Mf_time); % Convert MACA time to Matlab Datenum
        Mf_spd = hypot(Mf_u,Mf_v);
        wnddir_temp = (180/pi)*atan2(Mf_u,Mf_v);
        % Rotate winds to compass directions
        wnddir_temp = 90 - wnddir_temp;
        wnddir_temp(wnddir_temp<0)=wnddir_temp(wnddir_temp<0)+360;
        % switch to conventional coming from dir
        wnddir_temp = wnddir_temp+180;
        Mf_dir = wrap2360(wnddir_temp);
        clear wnddir_temp
end

% Concatenate Matrices 
M.u = cat(3,Mh_u,Mf_u);
clear Mh_u Mf_u
M.v = cat(3,Mh_v,Mf_v);
clear Mh_v Mf_v
M.time = cat(1,Mh_time,Mf_time);
clear Mh_time Mf_time
[M.lat,M.lon] = meshgrid(Mh_lat, Mh_lon);
clear Mh_lat Mh_lon 
M.spd = cat(3,Mh_spd,Mf_spd);
clear Mh_spd Mf_spd
M.dir = cat(3,Mh_dir,Mf_dir);
clear Mh_dir Mf_dir scenario


% Load OBS ----------------------------------------------------
fol_loc = 'C:\Users\ahooshmand\Desktop\Data_Forcings\station_data\gap_hourly'; % location of data
stations = {'bham_airport','mcchord_afb','vic_air','whidbey_nas'}; % station names
for m = 1:length(stations)
    file_load = strcat(fol_loc,'\',stations{m},'_hourly.mat');
    load(file_load);
    % Get Lat Lon values from different file
    fol_loc = 'C:\Users\ahooshmand\Desktop\Data_Forcings\station_data\nan_hourly';
    file_load = strcat(fol_loc,'\',stations{m},'_hourly.mat');
    temp_obs = load(file_load);
    % Grab meteo data for window of time
    obs(m).name = stations{m};
    obs(m).wndspd = wndspd;
    obs(m).wnddir = wnddir';
    obs(m).slp = slp;
    obs(m).time = time;
    obs(m).lat = temp_obs.lat;
    obs(m).lon = temp_obs.lon;
    clear airtemp slp time wnddir wndspd
end
load('WA_coast.mat');
% Find all shapes above threshold
wa_lat = [];
wa_lon = [];
thresh = 100 ;
for j = 1:length(wa_coast)
    temp_x = wa_coast(j).X;
    temp_y = wa_coast(j).Y;
    if length(temp_x) >= thresh %&& j ~= 3348 % 3348 is oregon
        for m = 1:length(temp_x);
            wa_lat(end+1) = temp_y(m);
            wa_lon(end+1) = temp_x(m);
        end
    end
end
clear j m temp_x temp_y thresh wa_coast
%% Find Closest MACA point to obs

for x = 1:length(obs)
    dist_list = (M.lon-obs(x).lon).^2 + (M.lat-obs(x).lat).^2;
    [~,I] = min(dist_list(:));
    temp_x = M.lon(I);
    temp_y = M.lat(I);
    [r(x),c(x)] = find(M.lon == temp_x & M.lat == temp_y);
end

%% Plot WA coastline with MACA points and OBS over top
clf 
p1 = plot(wa_lon, wa_lat);
hold on 
p2 = plot(M.lon, M.lat, 'r.');
legend([p2],'MACA Grid')
xlabel('Degrees Longitude')
ylabel('Degrees Latitude')
axis equal
hold on 
for x = 1:length(obs)
    plot(obs(x).lon, obs(x).lat,'ko','MarkerFaceColor','k')
    hold on 
%     plot(M.lon(r(x),c(x)),M.lat(r(x),c(x)),'r*')
%     hold on
end
xlim([-124, -122.4])
ylim([47.0 48.9])
%legend('MACA Points','Obs','Location','NorthWest')l
printFig(gcf, 'MACA_Position_With_Obs', [14 14], 'png', 150)

%% Plot Distributions for the stations with MACA 

%nbins = 15;
x_locs = 0:1:25;
% xbins = vecotrize
% Reshape obs to be multi-row x 1-column matrix of daily averages - Note
% column number is the station number from obs list
for s = 1:length(obs)
    temp_spd = obs(s).wndspd;
    day_vec = datenum(year(obs(s).time(1)),month(obs(s).time(1)),...
        day(obs(s).time(1))):1:datenum(year(obs(s).time(end)),...
        month(obs(s).time(end)),day(obs(s).time(end)));
    temp_avg = zeros(length(day_vec),1); % Initialize
    
    count = 0;
    for x = 1:round(length(temp_spd)/24)
        if x == 1
            i_st = 1;
            i_end = 24;
        else
            i_st = i_end +1;
            i_end = 24*x;
        end
        i_use = i_st:i_end;
        temp_avg(x) = nan_mean(temp_spd(i_use));
    end
    
    % Find where MACA is same as obs
    ma_inds = findnearest(M.time,day_vec(1)):findnearest(M.time,day_vec(end));
 
    % plot Distributions and save
    clf 
    histogram(temp_avg,x_locs,'Normalization','probability')
    hold on
    histogram(M.spd(:,:,ma_inds),x_locs,'Normalization','probability','FaceAlpha',0.4)
    legend('Obs - Daily Avg','MACA')
    xlabel('Wind Speed [m/s]')
    pause 
%     ylabel('Probability')
%     outname = sprintf('MACA_vs_Obs_Hist_%s',obs(s).name);
%     printFig(gcf, outname, [14 14], 'png', 150)
end
