%% Load in variables of interest
addpath C:\Functions_Matlab
clear all
close all
clc
% path = 'C:\Users\ahooshmand\Desktop\Xbeach\OwenBeach\ModelRuns\test\';
path = 'C:\Users\ahooshmand\Desktop\Xbeach\OwenBeach\ModelRuns\OB1_30points\';
filename = 'xboutput.nc';

xb_info = ncinfo([path filename]);

% Load Grid Info
X.x = load([path 'OB1_Xbeach_x.grd']);
X.y = load([path 'OB1_Xbeach_y.grd']);
X.z = load([path 'OB1_Xbeach_z.grd']);

X.point_zs = squeeze(ncread([path filename],'point_zs')); % Runup
X.point_xz = squeeze(ncread([path filename],'point_xz'));
X.point_yz = squeeze(ncread([path filename],'point_yz'));
X.zs = squeeze(ncread([path filename],'zs')); % Water Level

X.Hs = 4*sqrt(var(X.zs(:,1:end-1)'));
X.s = sqrt(X.x.^2+X.x.^2);
X.s = X.s - min(X.s);
%% Calculate Runup - Stockdon
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
swl = line([0 5000],[X.point_zs(1) X.point_zs(1)],'Color','k');
sup = line([0 5000],[setup setup],'LineStyle','--','Color','k');
xlabel('Time [quarter seconds]')
ylabel('Runup Height [m]')
set(gca,'fontsize',12)
legend([swl sup],'Still Water Level','Set up')
% legend([swl sup],'Still Water Level','Average Wave Setup','Location','Northwest','FontSize',14)
% printFig(gcf,'R2_Timeseries',[11 8.5],'png',300)

%% Plot cdf and get R2 from it 
clf
cdf = sort(pks,'ascend');
cdfy = linspace(0,1,length(cdf));
plot(cdf,cdfy)
grid on
hold on
l1 = plot([0 2.5],[.98 .98],'LineWidth',1.5);
l2 = line([1 1], [0 1],'Color','k','LIneWidth',1.2,'LineStyle','--');
l3 = line([setup setup],[0 1],'Color','m','LineWidth',1.2,'LineStyle','--');
box on
xlabel('Runup Elevation [m]','FontSize',14)
ylabel('Probability','FontSize',14)
set(gca,'fontsize',14)
ind = findnearest(0.98,cdfy);
legend([l1 l2 l3],'2% Exceedence','SWL','Average Wave Setup','Location','Northwest')
% printFig(gcf,'R2_CDF',[11 11],'png',300)

%% Calculate R2 and set up components from stockdon 

tides = load('E:\Abbas\Model_Met_Forcings\COOPS_tides\tacoma\Tacoma_Reconstruct_NAVD88');
tide_std=std(tides.twl);
upper_limit = X.point_zs

% Slope of shoreline
xx = [29.18 67.12];
yy = [2.455 -2.706];
beta = getSlope(xx,yy);
% beta = .28;

% Wave Forcing Statistics 
fo = .25;
Ho = 1.5;
k = getk(fo,Ho); 
Lo = 2*pi/k;

% Calculate R2 from cdf 
R2 = cdf(ind);

% Calculate R2 from stockdon
Stock_setup = 0.35*beta*sqrt(Ho*Lo);
Stock_swash = sqrt(Ho*Lo*(0.56*beta^2 + 0.004))/2;
Stock_R2 = 1.1*(Stock_setup + Stock_swash);

%% Plot Results 
clf
% Plot Land
fill([X.s X.s(end) X.s(end)], [X.z X.z(end) X.z(1)],[1 .9 .4],'LineStyle','none'); % Plot polygon of land
hold on 

% Add Model Results
wl = line([0 120], [X.point_zs(1) X.point_zs(1)],'LineWidth',2,'Color','b'); % Add water level 
mod_setup = line([0 120], [setup setup],'LineWidth',1,'Color','k'); % Add Model wave set up 
mod_r2 = line([0 120], [R2 R2],'LineWidth',2,'Color','k','LineStyle',':'); % Add Model R2

% Add Stockdon Calculations
s_setup = line([0 120], [X.point_zs(1)+Stock_setup X.point_zs(1)+Stock_setup],'LineWidth',1,'Color','r'); % Stockdon Setup 
s_r2 = line([0 120], [X.point_zs(1)+Stock_R2 X.point_zs(1)+Stock_R2],'LineWidth',1.5,'Color','r','LineStyle','--'); % Stockdon R2

% Aesthetics 
xlim([0 40])
ylim([-2 6])
xlabel('Cross Shore [m]')
ylabel('Elevation [m]')
set(gca,'FontSize',14)
legend([wl mod_setup mod_r2 s_setup s_r2],'SWL','Modeled Wave Setup','Modeled R2','Stockdon Wave Setup','Stockdon R2','Location','Northwest')

% printFig(gcf,'R2Comparison_ModVsStock',[11 8.5],'png',300)

%% Perform Same Analysis for Tillamook study 

% First calculate DWL - dynamic water level and Hmo
Dlow = -0.97;
n_s = .35*beta*sqrt(Ho*Lo); % Setup from stockdon
n_r = 0.06*sqrt(Ho*Lo); % IG swash from stockdon
DWL = X.point_zs(1) + 1.1 * (n_s + (n_r/2)) - Dlow; 
Hmo = DWL * 0.78;

% Calculate Iribarren Number
IB = beta/sqrt(Ho/Lo);

% Calculate Run-up
TawR2 = Hmo*(1.75*0.55*IB);