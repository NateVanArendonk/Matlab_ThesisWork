% THis script plots the WA nnrp point data with land mask and then the
% individual extracted points to show which ones were previously chosen 

load('WA_NNRP')

% Grab land and water points
land_lat = lat_m(land_mask);
land_lon = lon_m(land_mask);
water_lat = lat_m(~land_mask);
water_lon = lon_m(~land_mask);
return
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

% NNRP Points
np = dir('E:\Abbas\PS_COSMOS\Thesis_Modeling\Quantile_Correction\NNRP_PointData\*.mat');
for f = 1:length(np)
    f_open = strcat('E:\Abbas\PS_COSMOS\Thesis_Modeling\Quantile_Correction\NNRP_PointData\',np(f).name);
    n = load(f_open);
%     nn(f).lat = n.Position(2);
%     nn(f).lon = n.Position(1);
    nn(f).lat = n.lat;
    nn(f).lon = n.lon(1);
end
return


%% Plot 
clf
plot(wa_lon,wa_lat)
hold on 
p1 = plot(land_lon, land_lat, 'k.');
% Plot water Values
p2 = plot(water_lon, water_lat, 'r*');
xlabel('Degrees Longitude')
ylabel('Degrees Latitude')
axis equal
xlim([-123.3 -122.2])
ylim([46.9 48.4])
legend([p1 p2],'Land Points','Water Points','Location','NorthWest')
% printFig(gcf,'NNRPSouthPS_WaterLandPoints',[8.5 11],'png',300)
%% Plot all of PS Land Mask
clf
plot(wa_lon,wa_lat)
hold on 
p1 = plot(land_lon, land_lat, 'k.');
% Plot water Values
p2 = plot(water_lon, water_lat, 'r*');
xlabel('Degrees Longitude')
ylabel('Degrees Latitude')
axis equal
xlim([-124.75 -122.2])
ylim([46.9 49])
for ii = 1:length(nn)
    plot(nn(ii).lon,nn(ii).lat,'ob','MarkerSize',9)
end
legend([p1 p2],'Land Points','Water Points','Location','NorthWest')
printFig(gcf,'NNRP_ExtractedPoints',[8.5 11],'png',300)

