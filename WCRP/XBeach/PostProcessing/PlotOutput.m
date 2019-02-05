%% Load in variables of interest
addpath C:\Functions_Matlab
clear all
close all
clc

% path = 'C:\Users\ahooshmand\Desktop\Xbeach\OwenBeach\ModelRuns\test\';
path = 'C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\Tacoma\Owen\XBeach\OB_Runs\manning\50cmWaves\man25_contemp\';
grid_path = '../Owen_Beach_XbeachGrids/';
filename = 'xboutput.nc';

xb_info = ncinfo([path filename]);

% Load Grid Info
X.x = load([grid_path 'OB1_Xbeach_x.grd']);
X.y = load([grid_path 'OB1_Xbeach_y.grd']);
X.z = load([grid_path 'OB1_Xbeach_z.grd']);

X.point_zs = squeeze(ncread([path filename],'point_zs')); % Runup
X.point_xz = squeeze(ncread([path filename],'point_xz'));
X.point_yz = squeeze(ncread([path filename],'point_yz'));
X.zs = squeeze(ncread([path filename],'zs')); % Water Level

X.Hs = 4*sqrt(var(X.zs(:,1:end-1)')); % Gets rid of anomoly at end 
rho = 1027; %[kg/m3]
X.E = ((X.Hs.^2)*rho*9.81)/8; % Correct?

% Make cross shore 'S' transect
X.s = sqrt(X.x.^2+X.x.^2);
X.s = X.s - min(X.s);

% Make Cross shore 'S' transect for point Runup gauge data - BROKEN
X.point_s = sqrt(X.point_xz.^2+X.point_yz.^2);
X.point_s = X.point_s - min(X.point_s);
return
%% Plot Wave height with depth
clf
plot(X.s,X.Hs+X.zs(1,1))
hold on
plot(X.s,X.z)
grid on 
xlabel('Cross Shore [m]')
ylabel('Elevation [m]')
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
clf
% Calculate Set up 
time = 1:1:length(X.point_zs);
setup = nan_mean(X.zs);
setup = mean(setup);

% Find Peaks in Run up and plot 
[pks, locs] = findpeaks(X.point_zs,time,'MinPeakDistance',12,'MinPeakWidth',1);
plot(time,X.point_zs)
hold on
plot(locs,pks,'r*')
swl = line([0 length(time)],[X.point_zs(1) X.point_zs(1)],'Color','k');
sup = line([0 length(time)],[setup setup],'LineStyle','--','Color','k');
xlabel('Time [quarter seconds]')
ylabel('Runup Height [m]')
set(gca,'fontsize',12)
xlim([0 length(time)])
% legend([swl sup],'Still Water Level','Average Wave Setup','Location','Northwest','FontSize',14)
% printFig(gcf,'R2_Timeseries',[11 8.5],'png',300)

%% Plot cdf and get R2 from it 
clf
cdf = sort(pks,'ascend');
cdfy = linspace(0,1,length(cdf));
plot(cdf,cdfy)
grid on
hold on
l1 = plot([0 4],[.98 .98],'LineWidth',1.5);
l2 = line([X.zs(1,1) X.zs(1,1)], [0 1],'Color','k','LIneWidth',1.2,'LineStyle','--');
l3 = line([setup setup],[0 1],'Color','m','LineWidth',1.2,'LineStyle','--');
box on
xlabel('Runup Elevation [m]','FontSize',14)
ylabel('Probability','FontSize',14)
set(gca,'fontsize',14)
ind = findnearest(0.98,cdfy);
legend([l1 l2 l3],'2% Exceedence','SWL','Average Wave Setup','Location','Northwest')
% printFig(gcf,'R2_CDF',[11 11],'png',300)

% Calculate R2
R2 = cdf(ind);
%% Plot Results 
clf
% Plot Land
fill([X.s X.s(end) X.s(end)], [X.z X.z(end) X.z(1)],[1 .9 .4],'LineStyle','none'); % Plot polygon of land
hold on 

% Add Model Results
wl = line([0 200], [X.point_zs(1) X.point_zs(1)],'LineWidth',2,'Color','b'); % Add water level 
mod_setup = line([0 200], [setup setup],'LineWidth',1,'Color','k'); % Add Model wave set up 
mod_r2 = line([0 200], [R2 R2],'LineWidth',2,'Color','k','LineStyle',':'); % Add Model R2


% Aesthetics 
xlim([60 120])
ylim([-.5 6])
xlabel('Cross Shore [m]')
ylabel('Elevation [m]')
set(gca,'FontSize',14)
legend([wl mod_setup mod_r2 s_setup s_r2],'SWL','Modeled Wave Setup','Modeled R2','Stockdon Wave Setup','Stockdon R2','Location','Northwest')

printFig(gcf,'R2Comparison_ModVsStock_100grid',[11 8.5],'png',300)

