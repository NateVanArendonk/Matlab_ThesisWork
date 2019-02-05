% THis script plots the WA nnrp point data with land mask and then the
% individual extracted points to show which ones were previously chosen 

load('WA_NNRP')

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
return
MM = kml2struct('C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\Quantile_Correction\JDFMask.kml');
in = inpolygon(lon_m,lat_m,MM.Lon,MM.Lat);
jdf_x = lon_m(in);
jdf_y = lat_m(in);
% NNRP Points
np = dir('temp\*.mat');
for f = 1:length(np)
    f_open = strcat('temp\',np(f).name);
    n = load(f_open);
%     nn(f).lat = n.Position(2);
%     nn(f).lon = n.Position(1);
    nn(f).lat = n.lat;
    nn(f).lon = n.lon(1);
end

% Load Basin Masks
masks = dir('E:\Abbas\PS_COSMOS\Thesis_Modeling\KML\Basin_Masks\*.kml');
for f = 1:length(masks)
    Ma(f) = kml2struct(['E:\Abbas\PS_COSMOS\Thesis_Modeling\KML\Basin_Masks\' masks(f).name]);
end
return
%% Plot South PS Land Mask
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
legend([p1 p2],'Land Points','Water Points','Location','NorthWest')
printFig(gcf,'NNRPPugetSound_WaterLandPoints',[8.5 11],'png',300)

%% Plot NNRP and chosen points 
clf
plot(wa_lon,wa_lat,'Color',[.7 .7 .7])
hold on 
p1 = plot(land_lon, land_lat, 'k.');
% Plot water Values
p2 = plot(water_lon, water_lat, 'r*');
for ii = 1:length(nn)
    plot(nn(ii).lon,nn(ii).lat,'ob','MarkerSize',9)
end
xlabel('Degrees Longitude')
ylabel('Degrees Latitude')
axis equal
xlim([-124.75 -122.2])
ylim([46.9 49])
legend([p1 p2],'Land Points','Water Points','Location','NorthWest')
% printFig(gcf,'NNRP_ExtractedPoints',[11 11],'png',300)
%% Plot Land Masks on top of NNRP Points 
clf
plot(wa_lon,wa_lat,'Color',[.7 .7 .7])
xlim([-124.75 -122.2])
ylim([46.9 49])
hold on 
p1 = plot(land_lon, land_lat, 'k.');
% Plot water Values
p2 = plot(water_lon, water_lat, 'r*');
for ii = 1:length(nn)
    plot(nn(ii).lon,nn(ii).lat,'ob','MarkerSize',9)
end
% Add Land Masks
for ii = 1:length(Ma)
    pgon = polyshape(Ma(ii).Lon,Ma(ii).Lat);
    plot(pgon)
    alpha(0.3)
%     fprintf('Basin: %s',Ma(ii).Name)
%     pause
end
xlabel('Degrees Longitude')
ylabel('Degrees Latitude')
axis equal
xlim([-124.75 -122.2])
ylim([46.9 49])
legend([p1 p2],'Land Points','Water Points','Location','NorthWest')
% printFig(gcf,'NNRP_ExtractedPoints_WithBasins',[11 11],'png',300)

%% Confirm that basins and points sync up 

% Make key-value pairs for nnrp_points and polygons 
% keys = {'Anacortes.kml','BellinghamBay.kml','Blaine.kml','HoodCanal.kml',...
%     'LummiBay.kml','LummiIsland.kml','NorthCentral.kml','Olympia.kml','OlympiaNorth.kml','OlympiaWest.kml',...
%     'PadillaBay.kml','PortTownsend_Hood.kml','Skagit.kml','SouthCentral.kml',...
%     'Tacoma.kml','Whidbey_JDF.kml','WhidbeyAnacortes.kml','WhidbeyBasin.kml'...
%     'JDF1.kml','JDF2.kml','JDF3.kml','JDF4.kml','JDF5.kml','JDF6.kml','JDF7.kml','JDF8.kml','JDF9.kml','JDF10.kml','JDF11.kml'};
% values = {'Cypress_Guemes_NNRP_Point.mat','Chuckanut_NNRP_Point.mat',...
%     'Blaine_NNRP_Point.mat','Nor_Hood_Canal_NNRP_Point.mat','Lummi_Bay_NNRP_Point.mat',...
%     'Lummi_Island_NNRP_Point.mat','Kingston_NNRP_Point.mat','South_Tacoma_NNRP_Point.mat',...
%     'Nor_Olympia_NNRP_Point.mat','Middle_Olympia_NNRP_Point.mat','Padilla_Bay_NNRP_Point.mat',...
%     'Townsend_Whidbey_NNRP_Point.mat','Skagit_Bay_NNRP_Point.mat','blakeisland_corwith_westpoint.mat',...
%     'Tacoma_NNRP_Point.mat','JDF_Whidbey_NNRP_Point.mat','Deception_Pass_NNRP_Point.mat','Whidbey_Camano_NNRP_Point.mat'...
%     'JDF1_NNRP_Point.mat','Neah_Bay_NNRP_Point.mat','JDF3_NNRP_Point.mat','JDF4_NNRP_Point.mat','JDF5_NNRP_Point.mat'...
%     'JDF6_NNRP_Point.mat','JDF7_NNRP_Point.mat','JDF8_NNRP_Point.mat','JDF9_NNRP_Point.mat','Dungeness_NNRP_Point.mat','JDF11_NNRP_Point.mat'};
% ----------------- Key Value Pairs for Basin and NNRP --------------------
keys = {'Anacortes.kml','BellinghamBay.kml','Blaine.kml','Hood_Hook.kml',...
    'Hood_Kitsap.kml','Middle_Hood.kml','Port_Ludlow.kml','Nor_Hood.kml',...    
    'LummiBay.kml','LummiIsland.kml','NorthCentral.kml','Olympia.kml','OlympiaNorth.kml','OlympiaWest.kml',...
    'PadillaBay.kml','PortTownsend_Hood.kml','Skagit.kml','SouthCentral.kml',...
    'Tacoma.kml','Whidbey_JDF.kml','WhidbeyAnacortes.kml','WhidbeyBasin.kml'...
    'JDF1.kml','JDF2.kml','JDF3.kml','JDF4.kml','JDF5.kml','JDF6.kml','JDF7.kml',...
    'JDF8.kml','JDF9.kml','JDF10.kml','JDF11.kml','GuemesCypress_Mask1.kml',...
    'CypressSinclairGuemes_Mask1.kml','Guemes_Mask1.kml','LopezIsland_Mask1.kml',...
    'Orcas_Mask1.kml','SanJuan_Mask1.kml','SanJuan_Mask2.kml','SanJuan_mask3.kml',...
    'Stuart_Mask1.kml','Waldron_Mask1.kml','Waldron_Mask2.kml','Sucia_Mask1.kml',...
    'BoundaryBay_Mask.kml','PointRoberts_Mask.kml','Tsawwassen_Mask.kml',...
    'RichmondBC_Mask.kml','IonaBeachBC_Mask.kml','BurrardInlet_Mask.kml',...
    'BowenIslandBC_Mask.kml'};

values = {'Cypress_Guemes_NNRP_Point.mat','Chuckanut_NNRP_Point.mat',...
    'Blaine_NNRP_Point.mat','Hood_Hook_CORRECTED_NNRP_Point.mat',...
    'Hood_Kitsap_CORRECTED_NNRP_Point.mat','Middle_Hood_Canal_CORRECTED_NNRP_Point.mat',...
    'Port_Ludlow_CORRECTED_NNRP_Point.mat','Nor_Hood_Canal_CORRECTED_NNRP_Point.mat',...
    'Lummi_Bay_NNRP_Point.mat','Lummi_Island_NNRP_Point.mat','Kingston_NNRP_Point.mat',...
    'McNeil_Island_NNRP_Point.mat','Nor_Olympia_CORRECTED_NNRP_Point.mat',...
    'Middle_Olympia_CORRECTED_NNRP_Point.mat','Padilla_Bay_NNRP_Point.mat',...
    'Townsend_Whidbey_NNRP_Point.mat','Skagit_Bay_NNRP_Point.mat','Blake_Island_NNRP_Point.mat',...
    'Tacoma_CORRECTED_NNRP_Point.mat','JDF_Whidbey_NNRP_Point.mat','Deception_Pass_NNRP_Point.mat',...
    'Whidbey_Camano_NNRP_Point.mat','JDF1_NNRP_Point.mat','Neah_Bay_NNRP_Point.mat',...
    'JDF3_NNRP_Point.mat','JDF4_NNRP_Point.mat','JDF5_NNRP_Point.mat','JDF6_NNRP_Point.mat',...
    'JDF7_NNRP_Point.mat','JDF8_NNRP_Point.mat','JDF9_NNRP_Point.mat','Dungeness_NNRP_Point.mat',...
    'JDF11_NNRP_Point.mat','Cypress_Guemes_NNRP_Point.mat','Lummi_Island_NNRP_Point.mat',...
    'Padilla_Bay_NNRP_Point.mat','Lopez_South_NNRP_Point.mat','SJ_Islands_NNRP_Point.mat',...
    'SanJuan_South_NNRP_Point.mat','SanJuan_Island_NNRP_Point.mat','Stuart_Island_NNRP_Point.mat',...
    'SaturnaIsland_NNRP_Point.mat','Orcas_Island_NNRP_Point.mat','Saturna_Waldron_NNRP_Point.mat',...
    'Sucia_Island_NNRP_Point.mat','Boundary_Bay_NNRP_Point.mat','Point_Roberts_NNRP_Point.mat',...
    'Tsawwassen_NNRP_Point.mat','RichmondBC_NNRP_Point.mat','IonaBeach_BC_NNRP_Point.mat',...
    'BurrardInlet_BC_NNRP_Point.mat','BowenIsland_BC_NNRP_Point.mat'};
M = containers.Map(keys,values); % Creates a dictionary of key value pairs 

b = parula(length(unique(values)));
clf
plot(wa_lon,wa_lat,'Color',[.7 .7 .7])
axis equal
xlim([-125 -122.2])
ylim([46.9 49.5])
hold on 

colorNum = randi(45,length(b),1);
nameColors = cell(2,length(b));
nameColors(1,:) = unique(values);
for ii = 1:length(b)
    nameColors{2,ii} = colorNum(ii);
end


kml_fol = 'E:\Abbas\PS_COSMOS\Thesis_Modeling\KML\Basin_Masks\';
point_fol = 'E:\Abbas\PS_COSMOS\Thesis_Modeling\Quantile_Correction\NNRP_WaterPointData\';
for ii = 1:length(keys)
    K = kml2struct([kml_fol keys{ii}]);
    P = load([point_fol M(keys{ii})]);
    % Find what color of name should be and plot 
    name = M(keys{ii});
    I = strcmp(name,nameColors(1,:));
    I = find(I == 1);
    pgon = polyshape(K.Lon,K.Lat);
    plot(pgon,'FaceColor',b(nameColors{2,I},:),'FaceAlpha',0.8)
    plot(P.lon,P.lat,'o','Color',b(nameColors{2,I},:),'MarkerFaceColor',b(nameColors{2,I},:),'MarkerEdge','k') 
end
xlabel('Degrees Longitude')
ylabel('Degrees Latitude')
printFig(gcf,'BasinsWithNNRP_Point',[11 8.5],'png',300)

%% Color just the shoreline 

M = containers.Map(keys,values); % Creates a dictionary of key value pairs 

b = parula(length(unique(values)));
clf
plot(wa_lon,wa_lat,'Color',[.7 .7 .7])
axis equal
xlim([-125 -122.2])
ylim([46.9 49.5])
hold on 

colorNum = randi(45,length(b),1);
nameColors = cell(2,length(b));
nameColors(1,:) = unique(values);
for ii = 1:length(b)
    nameColors{2,ii} = colorNum(ii);
end

tic
kml_fol = 'E:\Abbas\PS_COSMOS\Thesis_Modeling\KML\Basin_Masks\';
point_fol = 'E:\Abbas\PS_COSMOS\Thesis_Modeling\Quantile_Correction\NNRP_WaterPointData\';
for ii = 1:length(keys)
    K = kml2struct([kml_fol keys{ii}]);
    % Find what color of name should be and plot 
    name = M(keys{ii});
    I = strcmp(name,nameColors(1,:));
    I = find(I == 1);
%     pgon = polyshape(K.Lon,K.Lat);
    in = inpolygon(wa_lon,wa_lat,K.Lon,K.Lat);
    plot(wa_lon(in),wa_lat(in),'.','Color',b(nameColors{2,I},:))
end
for ii = 1:length(keys)
    P = load([point_fol M(keys{ii})]);
    name = M(keys{ii});
    I = strcmp(name,nameColors(1,:));
    I = find(I == 1);
    plot(P.lon,P.lat,'o','Color',b(nameColors{2,I},:),'MarkerFaceColor',b(nameColors{2,I},:),'MarkerEdge','k','MarkerSize',10)
    hold on 
end
xlabel('Degrees Longitude')
ylabel('Degrees Latitude')
toc
printFig(gcf,'WashingtonShoreline_NNRP_PointCOLORED',[11 8.5],'png',300)

