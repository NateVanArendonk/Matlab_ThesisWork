%% Load Data
clearvars
% Add OET !!!!!!!!!!!!!
% Output Folder and File Location and Name
% ------- Will require Changing --------------------
nc_fol = 'E:\DelftFM\AndrewRuns\RoughnessTesting\Manning\r_316\';
nc_file = 'ps_2d_map.nc';
basemap = 'C:\Users\ahooshmand\Desktop\FM_Plotting\tiffs\bbay_1.mat';

% Get time info from MDU file
ref_date = datenum(2017,01,01); % Simulation Reference Date
sim_start = ref_date+(19440/(60*24)); % Start of simulation w.r.t ref_date
sim_end = ref_date+(129600/(60*24)); % End of simulation w.r.t ref_date
map_delay = 1166400/(60*60*24); % Dealy in writing map output from .mdu
map_start = sim_start+map_delay; % datenum start of map output
map_time = map_start:1/24:map_start+length(time)-1; % datenum of map output time


% ~ 5 seconds to load
G = dflowfm.readNet([nc_fol nc_file]); % Load in Grid Structure
D = dflowfm.readMap([nc_fol nc_file]); % Load in Data
bl = ncread([nc_fol nc_file],'FlowElem_bl',[1],[Inf]);
D.face.bl = bl;
% info = ncinfo([nc_fol nc_file]);
% T.datenum = nc_cf_time([nc_fol nc_file]); % Should load in time 

% load Individual Timeseries Data - 25 seconds to load
tic
vx = ncread([nc_fol nc_file],'ucx')'; % Velocities
vy = ncread([nc_fol nc_file],'ucy')';
time = ncread([nc_fol nc_file],'time');
wl = ncread([nc_fol nc_file],'waterdepth');
bl = ncread([nc_fol nc_file],'FlowElem_bl',[1],[Inf]);
toc
% Load Basemap
IM = load('C:\Users\ahooshmand\Desktop\FM_Plotting\tiffs\bbay_1UTM.mat');
% Note 4 base maps exist for bbay, bbay1,bbay2,bbay3,bbay4.  Higher the
% number means higher the resolution on the delta 
%% Plotting  in general
figure(1);clf 

% Default Plotting
dflowfm.plotMap(G,D) % Default plots WL 


%% Plot Bed Level
clf
imagesc(IM.xm,IM.ym,IM.im)
set(gca,'ydir','normal')
hold on 
dflowfm.plotMap(G,D,'parameter','bl') % Plot Bed Level
clrmap('C:\OET\matlab\applications\delft3d_matlab\colormaps\gist_earth.clrmap'); % Colormap that Deltares uses - feel free to change
caxis([-20 5]) % Bathy Limits

% xlim([-122.75 -122.4])
% ylim([48.55 48.85])
xlim([5.2*10^5 5.45*10^5])
ylim([5.377*10^6 5.412*10^6])
%% Plot WL with Velocity Vectors 
clf
imagesc(IM.xm,IM.ym,IM.im)
set(gca,'ydir','normal')
hold on 
dflowfm.plotMap(G,D,'parameter','dep')
colormap(parula)
caxis([-20 4])
hold on 
% quiver(G.face.FlowElem_x',G.face.FlowElem_y',D.face.u,D.face.v,'w')
for i = 1:10
    quiver(G.face.FlowElem_x',G.face.FlowElem_y',vx(i,:)',vy(i,:)','w')
    pause(0.25)
%     xlim([-122.75 -122.4])
%     ylim([48.55 48.85])
xlim([5.2*10^5 5.45*10^5])
ylim([5.377*10^6 5.412*10^6])
end
% xlim([-122.75 -122.4])
% ylim([48.55 48.85])
xlim([5.2*10^5 5.45*10^5])
ylim([5.377*10^6 5.412*10^6])
%% Plot WL through time 
% NOTE FOR SOME REASON WITH PLOTTING YOU CANNOT STOP THE LOOP WITH THIS OET
% PLOTTING SO TEST WITH SHORT TIME PERIODS!
clf
imagesc(IM.xm,IM.ym,IM.im)
set(gca,'ydir','normal')
hold on 
for i = 500:600
    dflowfm.plotMap(G,wl(:,i))
    colormap(parula)
    pause(0.25)
    caxis([-4 4])
%     xlim([-122.75 -122.4])
%     ylim([48.55 48.85])
xlim([5.2*10^5 5.45*10^5])
ylim([5.377*10^6 5.412*10^6])
end
% xlim([-122.75 -122.4]) % limits for bbay1
% ylim([48.55 48.85])
xlim([5.2*10^5 5.45*10^5])
ylim([5.377*10^6 5.412*10^6])
