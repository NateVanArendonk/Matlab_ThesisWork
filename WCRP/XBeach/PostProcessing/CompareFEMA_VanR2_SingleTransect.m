% This code will compare TWL estimates between FEMA, Site Workshop, Van der
% Meer and Xbeach 
clearvars
addpath C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\RunUp\RunupFormulas
% First load in the DEM 
E = load('E:\Abbas\Modeling Resources\PS_DEM\Ruston_Way\OwenBeach_CONED_Subset.mat');
load('C:\Users\ahooshmand\Desktop\PS_COSMOS\Salish_Model_Resources\WA_Spatial_Data\WA_coast_UTM');

% Subset even further to be owen beach 
inds = E.y_s<= 5.241*10^6 & E.x_s <= 5.366*10^5 & E.x_s >= 5.35*10^5;
E.x_s = E.x_s(inds);
E.y_s = E.y_s(inds);
E.z_s = E.z_s(inds);

% Load in Waves 
W = load('C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\Hindcast\Ruston_LUT\OwenRuston_wave_hindcast_SLR0.0FT.mat');

% Make a gridded interped surface of DEM
myX = E.x_s(1):1:E.x_s(end);
myY = E.y_s(end):1:E.y_s(1);
tic
[mx,my] = meshgrid(myX,myY);
mz = griddata(E.x_s,E.y_s,E.z_s,mx,my);
toc

% Load in FEMA FIRM 100 year 
F = load('E:\Abbas\WCRP\Tacoma\FEMA\Ruston_1Pct_Annual_Chance_Flood_Line');
% convert to UTM
[F.x,F.y] = deg2utm(F.lat,F.lon);
% Load in  Site Workshop study 
siteMHHW = 2.83; %[m] - NAVD
site2070 = 3.23; %[m] - NAVD - 50% scenario by 2070, 1.3 feet SLR 
site2120 = 3.69; %[m] - NAVD - 50% scenario by 2120, 2.8 feet SLR

% ---------- Location of KKL - Load and Convert to UTM 
kml_fol = 'C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\KML\RunupTransects\OwenBeach\';
kml_nm = 'OB1.kml';
temp = kml2struct([kml_fol kml_nm]); % load in KML of Transect
T.lat = temp.Lat;
T.lon = temp.Lon;
[T.x_utm,T.y_utm] = deg2utm(T.lat,T.lon); % Convert to utm
T.name = kml_nm;
clear kml_fol kmls kml_nm

% Find closest Fema point to end of transect 
dist_list = sqrt((T.x_utm(2) - F.x).^2 + (T.y_utm(2)-F.y).^2);
[~,I] = min(dist_list);
% Grab elevation at that point
[r,c] = findNearestGridPoint(mx,my,F.x(I),F.y(I));
F100 = mz(r,c);

% Find Closest Wave point to start of transect 
% dist_list = sqrt((T.x_utm(1) - W.x).^2 + (T.y_utm(1)-W.y).^2);
% [~,I] = min(dist_list);
% hs = W.hs_ts(I,:);
% tp = W.tp_ts(I,:);

% Load in Vandermeer Contemporary conditions
V = load('OB1_R2_rough80.mat');

twl = V.V2+W.twl';
return

%%  Plot FEMA and R2

clf
 % Load in Geotiff
IM = load('C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\GeoTiffs\OwenBeach\OB_utm_zoom2.mat');
imagesc(IM.xm,IM.ym,IM.im)
set(gca,'ydir','normal')
hold on

% Add Contours 
contour_levels = (0:.25:8);
[C,h] = contour(mx,my,mz,contour_levels); % use contourf to plot filled in bathy 
% c = colorbar;
% c.Label.String = 'Elevation [m - NAVD88]';
h.LineColor = [.7 .7 .7];
xlabel('Northing [m]');
ylabel('Easting [m]');

% Label Contours 
v = 0:1:8;
clabel(C,h,v)

% Add in DEM and make transparent 
% p = pcolor(mx,my,mz);
% shading interp
% set(p,'facealpha',0.3)
% caxis([0 8])


% Add Transect 
plot(T.x_utm,T.y_utm,'w')

xlim([5.3560*10^5, 5.357*10^5])
ylim([5.24001*10^6, 5.24008*10^6])


% Add Fema 
plot(F.x,F.y,'r','LineWidth',1.5)

% Add R2 contour 
cR2 = [max(twl) max(twl)];
cr2 = contour(mx,my,mz,cR2,'b');
cr2(:,1) = [];
temp = cr2; clear cr2
cr2.x = temp(1,:);
cr2.y = temp(2,:);
cr2 = plot(cr2.x(1:end-100),cr2.y(1:end-100),'b','LineWidth',1.5);
femaText = sprintf('FEMA: %.2f m NAVD88',F100);
vText = sprintf('TAW: %.2f m NAVD88',max(twl));
text(5.35675*10^5,5.24005*10^6,femaText,'FontSize',14,'Color','W')
text(5.35675*10^5,5.240047*10^6,vText,'FontSize',14,'Color','W')

set(gca,'FontSize',14)
% printFig(gcf,'OB1_FEMA_R2',[15 11],'png',300)


%% Plot of varying roughness 

% ---------------- Add R2 contour 
cR2 = [max(twl) max(twl)];
cr2 = contour(mx,my,mz,cR2,'b');
cr2(:,1) = [];
temp = cr2; clear cr2
cr2.x = temp(1,:);
cr2.y = temp(2,:);

% ------------ Add different roughness
cS1 = [max(twlR1),max(twlR1)];
cs1 = contour(mx,my,mz,cS1,'m');
cs1(:,1) = [];
temp = cs1; clear cs1
cs1.x = temp(1,:);
cs1.y = temp(2,:);

% ------------ Add other roughness
cS2 = [max(twlR2) max(twlR2)];
cs2 = contour(mx,my,mz,cS2,'m');
cs2(:,1) = [];
temp = cs2; clear cs2
cs2.x = temp(1,:);
cs2.y = temp(2,:);


clf
 % Load in Geotiff
IM = load('C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\GeoTiffs\OwenBeach\OB_utm_zoom2.mat');
imagesc(IM.xm,IM.ym,IM.im)
set(gca,'ydir','normal')
hold on

% Add Contours 
contour_levels = (0:.25:8);
[C,h] = contour(mx,my,mz,contour_levels); % use contourf to plot filled in bathy 
% c = colorbar;
% c.Label.String = 'Elevation [m - NAVD88]';
h.LineColor = [.7 .7 .7];
xlabel('Northing [m]');
ylabel('Easting [m]');

% Label Contours 
v = 0:1:8;
clabel(C,h,v)

% Add Transect 
plot(T.x_utm,T.y_utm,'w')

% Add Fema 
plot(F.x,F.y,'r','LineWidth',1.5)

% Add Runup values 
cr2 = plot(cr2.x(1:end-100),cr2.y(1:end-100),'b','LineWidth',1.5);
cs1 = plot(cs1.x(1:end-100),cs1.y(1:end-100),'m','LineWidth',1.5);
cs2 = plot(cs2.x(1:end-100),cs2.y(1:end-100),'c','LineWidth',1.5);


xlim([5.3560*10^5, 5.357*10^5])
ylim([5.24001*10^6, 5.24008*10^6])

%% Add SLR Scenarios for R2

% ---------------- Add R2 contour 
cR2 = [R2, R2];
cr2 = contour(mx,my,mz,cR2,'b');
cr2(:,1) = [];
temp = cr2; clear cr2
cr2.x = temp(1,:);
cr2.y = temp(2,:);

% ------------ Add SLR scenario #1: 50% by 2070 of 1.3' 
twl = V1.V2 + V1.twl+(.3048*1.3); % GET RID OF SLR VALUE ONCE YOU RERUN
[pks, locs] = findpeaks(twl,W.time);

cdf = sort(pks,'ascend');
cdfy = linspace(0,1,length(cdf));

% Find 2% run up elevation 
ind = findnearest(0.98,cdfy);
s1R2 = cdf(ind);
% Get contour of slr value 
cS1 = [s1R2, s1R2];
cs1 = contour(mx,my,mz,cS1,'m');
cs1(:,1) = [];
temp = cs1; clear cs1
cs1.x = temp(1,:);
cs1.y = temp(2,:);

% ------------ Add SLR scenario #2: 50% by 2012 of 2.8' 
twl = V2.V2 + V2.twl+(.3048*2.8); % GET RID OF SLR VALUE ONCE YOU RERUN
[pks, locs] = findpeaks(twl,W.time);

cdf = sort(pks,'ascend');
cdfy = linspace(0,1,length(cdf));

% Find 2% run up elevation 
ind = findnearest(0.98,cdfy);
s2R2 = cdf(ind);
% Get contour of slr value 
cS2 = [s2R2, s2R2];
cs2 = contour(mx,my,mz,cS2,'m');
cs2(:,1) = [];
temp = cs2; clear cs2
cs2.x = temp(1,:);
cs2.y = temp(2,:);


clf
 % Load in Geotiff
IM = load('C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\GeoTiffs\OwenBeach\OB_utm_zoom2.mat');
imagesc(IM.xm,IM.ym,IM.im)
set(gca,'ydir','normal')
hold on

% Add Contours 
contour_levels = (0:1:8);
[C,h] = contour(mx,my,mz,contour_levels); % use contourf to plot filled in bathy 
% c = colorbar;
% c.Label.String = 'Elevation [m - NAVD88]';
h.LineColor = [.7 .7 .7];
xlabel('Northing [m]');
ylabel('Easting [m]');

% Label Contours 
v = 0:1:8;
clabel(C,h,v)

% Add Transect 
plot(T.x_utm,T.y_utm,'w')

% Add Fema 
plot(F.x,F.y,'r','LineWidth',1.5)

% Add Runup values 
cr2 = plot(cr2.x(1:end-100),cr2.y(1:end-100),'b','LineWidth',1.5);
cs1 = plot(cs1.x(1:end-100),cs1.y(1:end-100),'m','LineWidth',1.5);
cs2 = plot(cs2.x(1:end-100),cs2.y(1:end-100),'c','LineWidth',1.5);


xlim([5.3560*10^5, 5.357*10^5])
ylim([5.24001*10^6, 5.24008*10^6])