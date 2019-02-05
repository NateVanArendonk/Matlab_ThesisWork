clearvars
%OBS
addpath C:\Functions_Matlab\time

%stations = dir('C:\Users\ahooshmand\Desktop\Data_Forcings\station_data\Longer_Record_Data\*.mat'); % Longer Record Data
stations = dir('C:\Users\ahooshmand\Desktop\Data_Forcings\station_data\nan_hourly\*.mat'); % All Stations
for f = 1:length(stations)
    f_open = strcat('C:\Users\ahooshmand\Desktop\Data_Forcings\station_data\nan_hourly\',stations(f).name);
    O = load(f_open);
    if length(O.lat) > 1
        obs(f).lat = O.lat(1);
        obs(f).lon = -1*O.lon(1);
    else
        obs(f).lat = O.lat;
        obs(f).lon = O.lon;
    end
    
    % Grab Coverage
    st_yr = year(O.time(1));
    ed_yr = year(O.time(end));
    obs(f).cov = abs(ed_yr - st_yr);
end

clear O
% NNRP Points
np = dir('NNRP_Points\*.mat');
for f = 1:length(np)
    f_open = strcat('NNRP_Points\',np(f).name);
    n = load(f_open);
    nn(f).lat = n.Position(2);
    nn(f).lon = n.Position(1);
end
fol_loc = 'C:\Users\ahooshmand\Desktop\Data_Forcings\Wind_Products\ExtractNNRP\';
lat_m = ncread(strcat(fol_loc,'geo_em.d02.nc'),'XLAT_M');
lon_m = ncread(strcat(fol_loc,'geo_em.d02.nc'),'XLONG_M');
land_mask = ncread(strcat(fol_loc,'geo_em.d02.nc'),'LANDMASK');
land_mask = logical(land_mask);
% Grab land and water points
land_lat = lat_m(land_mask);
land_lon = lon_m(land_mask);
water_lat = lat_m(~land_mask);
water_lon = lon_m(~land_mask);

% Load Coastline
load('WA_coast.mat');
% Find all shapes above threshold
wa_lat = [];
wa_lon = [];
thresh = 500 ;
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
%NN = load('NNRP_GRID_POINTS');
clear f j m thresh temp_x temp_y wa_coast f_open
% Plotting
clf
plot(wa_lon,wa_lat)
hold on 
% Plot obs
for x = 1:length(obs)
    coverage = obs(x).cov;
    disp(coverage)
    if coverage <= 15
        plot(obs(x).lon,obs(x).lat,'hk','MarkerFaceColor','k')
    elseif coverage >= 16 && coverage <= 24
        plot(obs(x).lon,obs(x).lat,'dk','MarkerFaceColor','k')
    elseif coverage >= 25
        plot(obs(x).lon,obs(x).lat,'ok','MarkerFaceColor','k')
    end
    hold on
end
% Plot Land Values
plot(land_lon, land_lat, 'k.')
% Plot water Values
plot(water_lon, water_lat, 'r*')


%plot(NN.lon,NN.lat,'.r')

% for x = 1:length(nn)
%     plot(nn(x).lon,nn(x).lat,'*r');
%     hold on
% end 
xlabel('Degrees Longitude')
ylabel('Degrees Latitude')
axis equal
xlim([-123.3 -122.2])
ylim([46.9 48.4])
% xlim([-124 -122.2])
% ylim([46.9 50])
return
%% Plot South PS Land Mask
clf
plot(wa_lon,wa_lat)
hold on 
p1 = plot(land_lon, land_lat, 'k.')
% Plot water Values
p2 = plot(water_lon, water_lat, 'r*')
xlabel('Degrees Longitude')
ylabel('Degrees Latitude')
axis equal
xlim([-123.3 -122.2])
ylim([46.9 48.4])
legend([p1 p2],'Land Points','Water Points','Location','NorthWest')
printFig(gcf,'NNRPSouthPS_WaterLandPoints',[8.5 11],'png',300)
%% Plot all of PS Land Mask
clf
plot(wa_lon,wa_lat)
hold on 
p1 = plot(land_lon, land_lat, 'k.')
% Plot water Values
p2 = plot(water_lon, water_lat, 'r*')
xlabel('Degrees Longitude')
ylabel('Degrees Latitude')
axis equal
xlim([-124.75 -122.2])
ylim([46.9 49])
legend([p1 p2],'Land Points','Water Points','Location','NorthWest')
printFig(gcf,'NNRPPugetSound_WaterLandPoints',[8.5 11],'png',300)