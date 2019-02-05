%% Load in variables of interest
addpath C:\Functions_Matlab
clear all
close all
clc

% path = 'C:\Users\ahooshmand\Desktop\Xbeach\OwenBeach\ModelRuns\test\';
path = '../OB_Runs/RoughnessTesting/manning/';
grid_path = '../Owen_Beach_XbeachGrids/';
filename = '/xboutput.nc';

xbr = dir('../OB_Runs/RoughnessTesting/manning/*');
xbr(1:2) = []; % Gets rid of 2 weird results 

xb_info = ncinfo([path 'r1' filename]); % Just get info from on eof the runs 

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
%% Plot waterlevel
% wave_setup = nanmean(zs);
figure(4)
for ii = 1:size(X.zs,1)
    % ------- Plot Profiles Only
    % plot(X.s,X.z,'-','Linewidth',3,'Color',[1 .9 .4]) % Plot initial profile - keeping just in case you want to get rid of fill
    % hold on
    % plot(X.s,X.zs(:,ii),'b-','Linewidth',2) % plot water level for time step
    
    % ------- Plot Filled Polygons
    fill([X.s X.s(end) X.s(1)], [X.zs(:,ii)' X.z(1) X.z(1)],'b','LineStyle','none') % Plot polygon of water
    hold on
    fill([X.s X.s(end) X.s(end)], [X.z X.z(end) X.z(1)],[1 .9 .4],'LineStyle','none') % Plot polygon of land
    grid on 
    hold off
    title(['t = ' sprintf('%04d',ii) ' seconds'])
    xlabel('Cross Shore [m]')
    ylabel('Wave height [m]')
    axis([min(X.s) max(X.s) -3 6])
    pause(0.25)
    % set(gcf, 'Position', [100 100 800 600]);
end

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

% Add Model Results
wl = line([0 200], [X(1).point_zs(1) X(1).point_zs(1)],'LineWidth',2,'Color','b','LineStyle','--'); % Add water level 
b = parula(9);
for ii = 1:length(X)
    mod_r2(ii) = line([0 200], [X(ii).R2 X(ii).R2],'LineWidth',2,'Color',b(ii,:),'LineStyle','-'); % Add Model R2
end

% Add Vandermeer R2 value 
v = line([0 200], [3.4594 3.4594],'LineWidth',2,'Color','r','LineStyle','--'); 

% Aesthetics 
xlim([89 99])
ylim([3.1 3.55])
xlabel('Cross Shore [m]')
ylabel('Elevation [m]')
set(gca,'FontSize',14)
grid on 

lgd = legend(mod_r2,arrayfun(@(x)(x.rough),X,'un',0),...
    'interpreter','none',...
    'location','northeast');

printFig(gcf,'ManningRoughnessTesting',[15 11],'png',300)
