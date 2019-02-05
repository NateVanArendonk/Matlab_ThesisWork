%% Load in variables of interest
addpath C:\Functions_Matlab
clear all
close all
clc

waves = 'small'; % Load in 50cm or 15cm wave case 
switch waves 
    case 'small'
        path = '../OB_Runs/manning/15cmWaves/';
        grid_path = '../Owen_Beach_XbeachGrids/';
        filename = '/xboutput.nc';
        
        xbr = dir('../OB_Runs/manning/15cmWaves/*');
        xbr(1:2) = []; % Gets rid of 2 weird results
        
        xb_info = ncinfo([path 'man25_SLR0.40' filename]); % Just get info from on eof the runs
    case 'large'
        path = '../OB_Runs/manning/50cmWaves/';
        grid_path = '../Owen_Beach_XbeachGrids/';
        filename = '/xboutput.nc';
        
        xbr = dir('../OB_Runs/manning/50cmWaves/*');
        xbr(1:2) = []; % Gets rid of 2 weird results
        
        xb_info = ncinfo([path 'man25_SLR0.40' filename]); % Just get info from on eof the runs
end

for ii = 1:length(xbr)
X(ii).rough = xbr(ii).name;
    
% Load Grid Info
X(ii).x = load([grid_path 'OB1_Xbeach_x.grd']);
X(ii).y = load([grid_path 'OB1_Xbeach_y.grd']);
X(ii).z = load([grid_path 'OB1_Xbeach_z.grd']);

% Load wave info
X(ii).point_zs = squeeze(ncread([path xbr(ii).name filename],'point_zs')); % Runup
X(ii).point_xz = squeeze(ncread([path xbr(ii).name filename],'point_xz'));
X(ii).point_yz = squeeze(ncread([path xbr(ii).name filename],'point_yz'));
X(ii).zs = squeeze(ncread([path xbr(ii).name filename],'zs')); % Water Level

X(ii).Hs = 4*sqrt(var(X(ii).zs(:,1:end-1)')); % Gets rid of anomoly at end 
rho = 1027; %[kg/m3]
X(ii).E = ((X(ii).Hs.^2)*rho*9.81)/8; % Correct?

% Make cross shore 'S' transect
X(ii).s = sqrt(X(ii).x.^2+X(ii).x.^2);
X(ii).s = X(ii).s - min(X(ii).s);

% Make Cross shore 'S' transect for point Runup gauge data - BROKEN
X(ii).point_s = sqrt(X(ii).point_xz.^2+X(ii).point_yz.^2);
X(ii).point_s = X(ii).point_s - min(X(ii).point_s);
end

% Get FEMA 100 yr levels
% First load in the DEM 
E = load('E:\Abbas\Modeling Resources\PS_DEM\Ruston_Way\OwenBeach_CONED_Subset.mat');
load('C:\Users\ahooshmand\Desktop\PS_COSMOS\Salish_Model_Resources\WA_Spatial_Data\WA_coast_UTM');

% Subset even further to be owen beach 
inds = E.y_s<= 5.241*10^6 & E.x_s <= 5.366*10^5 & E.x_s >= 5.35*10^5;
E.x_s = E.x_s(inds);
E.y_s = E.y_s(inds);
E.z_s = E.z_s(inds);
% Grid that B
myX = E.x_s(1):1:E.x_s(end);
myY = E.y_s(end):1:E.y_s(1);
tic
[mx,my] = meshgrid(myX,myY);
mz = griddata(E.x_s,E.y_s,E.z_s,mx,my);
toc
% ---------- Location of Transect - Load and Convert to UTM 
kml_fol = 'C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\KML\RunupTransects\OwenBeach\';
kml_nm = 'OB1.kml';
temp = kml2struct([kml_fol kml_nm]); % load in KML of Transect
T.lat = temp.Lat;
T.lon = temp.Lon;
[T.x_utm,T.y_utm] = deg2utm(T.lat,T.lon); % Convert to utm
T.name = kml_nm;
clear kml_fol kmls kml_nm

% Load in FEMA FIRM 100 year 
F = load('C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\Tacoma\FEMA\Ruston_1Pct_Annual_Chance_Flood_Line');
% convert to UTM
[F.x,F.y] = deg2utm(F.lat,F.lon);
% Find closest Fema point to end of transect 
dist_list = sqrt((T.x_utm(2) - F.x).^2 + (T.y_utm(2)-F.y).^2);
[~,I] = min(dist_list);
% Grab elevation at that point
[r,c] = findNearestGridPoint(mx,my,F.x(I),F.y(I));
F100 = mz(r,c);

return
%% Calculate Runup - Stockdon methodology 

% Calculate Set up 
time = 1:1:length(X(1).point_zs);

% Find Peaks in Run up and plot 
for ii = 1:length(X)
[pks, locs] = findpeaks(X(ii).point_zs,time,'MinPeakDistance',12,'MinPeakWidth',1);
cdf = sort(pks,'ascend');
cdfy = linspace(0,1,length(cdf));
ind = findnearest(0.98,cdfy);
% Calculate R2
X(ii).R2 = cdf(ind);
end
%% Plot Results 
clf
% Plot Land
fill([X(1).s X(1).s(end) X(1).s(end)], [X(1).z X(1).z(end) X(1).z(1)],[1 .9 .4],'LineStyle','none'); % Plot polygon of land
hold on 

b = parula(100);
switch waves 
    case 'small'
        v = line([0 200], [3.87 3.87],'LineWidth',2,'Color',b(1,:),'LineStyle','--');
        
        % Add contemporary Xbeach
        x1 = line([0 200],[X(3).R2 X(3).R2],'LineWidth',2,'Color','r','LineStyle','-');
        
        % 1.3 ft slr
        s1 =  line([0 200], [X(1).R2 X(1).R2],'LineWidth',2,'Color',b(20,:),'LineStyle','-'); % Add Model R2
        
        % Add FEMA 100 year line
        f = line([0 200], [F100 F100],'LineWidth',2,'Color',b(50,:),'LineStyle','--');
        
        % 2.8 ft slr
        s2 =  line([0 200], [X(2).R2 X(2).R2],'LineWidth',2,'Color',b(69,:),'LineStyle','-'); % Add Model R2
        xlim([50 100])
        ylim([2 5])
        lgd = legend([s2,s1,f,v,x1],'2.8ft SLR - XBeach','1.3ft SLR - XBeach','FEMA 100yr line','Vandermeer R2 - Contemporary','XBeach - Contemporary','Location','NorthEast');
        lgd.FontSize = 14;
        fname = 'Xbeach_SLR_Comparison_SmallWaveEvent';
    case 'large'
        v = line([0 200], [3.4594 3.4594],'LineWidth',2,'Color',b(1,:),'LineStyle','--');
        
        % Add contemporary Xbeach
        x1 = line([0 200],[X(3).R2 X(3).R2],'LineWidth',2,'Color','r','LineStyle','-');
        
        % 1.3 ft slr
        s1 =  line([0 200], [X(1).R2 X(1).R2],'LineWidth',2,'Color',b(20,:),'LineStyle','-'); % Add Model R2
        
        % Add FEMA 100 year line
        f = line([0 200], [F100 F100],'LineWidth',2,'Color',b(50,:),'LineStyle','--');
        
        % 2.8 ft slr
        s2 =  line([0 200], [X(2).R2 X(2).R2],'LineWidth',2,'Color',b(69,:),'LineStyle','-'); % Add Model R2
        xlim([50 100])
        ylim([2 5])
        lgd = legend([s2,f,s1,x1,v],'2.8ft SLR - XBeach','FEMA 100yr line','1.3ft SLR - XBeach','XBeach - Contemporary','Vandermeer R2 - Contemporary','Location','NorthEast');
        lgd.FontSize = 14;
        fname = 'Xbeach_SLR_Comparison_LargeWaveEvent';
end


xlabel('Cross Shore [m]')
ylabel('Elevation [m]')
set(gca,'FontSize',14)
grid on 



printFig(gcf,fname,[15 11],'png',300)

