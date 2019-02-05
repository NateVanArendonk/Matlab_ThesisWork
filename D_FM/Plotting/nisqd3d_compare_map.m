% nisqd3d_compare_map - Read and plot flow and sediment map-data of
% several runs

%% Load OET
cd('C:\Users\Stefan\Documents\toolbox\OET\matlab\');
run oetsettings.m
cd('C:\Users\Stefan\Documents\scripts\');

%%
clear; close all; clc

%% Initialization
% Define folders
fldr.matlab  = 'c:\Users\Stefan\Documents\scripts\';
fldr.d3d     = 'c:\Users\Stefan\Documents\delft3d\results\';
fldr.figures = 'c:\Users\Stefan\Documents\figures\';
fldr.nmdir   = 'sp_mtrx1\';
% Create fldr.model directories and fldr.model names
div = {'REF' 'A' 'B' 'C'};
q   = [43 110 290 462];
sl  = [0 1];
for i = 1:length(div)
    for j = 1:length(q)
        for k = 1:length(sl)
        mod{i,j,k}   = [div{j} '_q' sprintf('%d',q(i)) '_sed1_sl' sprintf('%d',sl(k)) '\'];
        modnm{i,j,k} = [div{j} '' sprintf('%d',find(q(i) == q)) '-sl' sprintf('%d',sl(k))];
        end
    end
end
fldr.model = mod(:);
nmmodel    = modnm(:);
% Read D3D trim-files and 1 trih-file (for cross-sections and stations)
for i = 1:length(fldr.model)
    trim(i) = qpfopen([fldr.d3d fldr.nmdir fldr.model{i} 'trim-nisq.dat']);
    trih = qpfopen([fldr.d3d fldr.nmdir fldr.model{1} 'trih-nisq.dat']);
end
% Read time series trih-file and trim-file
spinup          = 7.2917;    % Spin-up (d)
trih_timedat    = qpread(trih(1),'water level','times');
trih_tidx_first = find((trih_timedat-trih_timedat(1))>=spinup,1,'first');
trih_time       = qpread(trih(1),'water level','times',(trih_tidx_first-1):length(trih_timedat));
trim_time       = qpread(trim(1),'water level','times',0);
% Read mgrid coordinates: X and Y
mgrid = qpread(trim(1),0,'hydrodynamic grid','griddata');
X = (mgrid.X(:,:,1)')./1000;
Y = (mgrid.Y(:,:,1)')./1000;
% Read coordinates polygons: diversion and restoration
cd([fldr.d3d fldr.nmdir fldr.model{1}])
[Xp, Yp] = landboundary('read','diversion.pol'); % Polygon defined to prevent large vectors overshadowing the diversion area
[Xr, Yr] = landboundary('read','restoration2009.pol'); % Polygon marking the restoration area (polder area)
cd(fldr.matlab)
Xp = (Xp./1000)';
Yp = (Yp./1000)';
in = inpolygon(X,Y,Xp,Yp);
Xr = (Xr./1000)';
Yr = (Yr./1000)';
% Stations
nmstn    = station(trih(1),'name');
coordstn = station(trih(1),'coordXY')./1000;
% Cross-sections
nmcrs  = crosssec(trih(1),'names');
coordcrs = crosssec(trih(1),'coordXY')./1000;

if ~exist([fldr.figures fldr.nmdir 'Delta_VelocityMaps'],'dir');
    mkdir([fldr.figures fldr.nmdir 'Delta_VelocityMaps']);
end

if ~exist([fldr.figures fldr.nmdir 'Delta_TransportMaps'],'dir');
    mkdir([fldr.figures fldr.nmdir 'Delta_TransportMaps']);
end

clear mod modnm i j k

%% 1 Delta Mean Velocity

% Select model runs (never select more than 2 runs)
div_s = {'B','C'}; % Select from 'REF', 'A', 'B', or 'C'
q_s   = [110];       % Only select 1 discharge
sl_s  = [0];         % Select from '0' or '1'
% Adjust size scale vector and color-axis based on discharge class 
switch q_s
    case {43 110}
        u_vec = 0.25;   % Magnitude (m/s)
        v_vec = 0;
        scale = 0.2;   % Scaling
        colax = [-0.25 0.25];        
    case {290 462}
        u_vec = 0.5;   % Magnitude (m/s)
        v_vec = 0;
        scale = 0.1;   % Scaling
        colax = [-0.5 0.5];
end
% Domains
    % Delta
crdfig.x1 = 5.20E2;
crdfig.x2 = 5.24E2;
crdfig.y1 = 52.13E2;
crdfig.y2 = 52.17E2;
    % Diversion
crdfig.x1div = 5.216E2;
crdfig.x2div = 5.225E2;
crdfig.y1div = 52.138E2;
crdfig.y2div = 52.147E2;
% Select model runs
tridx = [];
for i = 1:length(div_s) 
    for j = 1:length(q_s)
        for k = 1:length(sl_s)
        tridx(i,j,k)  = find(strcmpi([div_s{i} '' sprintf('%d',find(q_s(j) == q)) '-sl' sprintf('%d',sl_s(k))],nmmodel));   
        end
    end
end
tridx = sort(tridx(:));
names = nmmodel(tridx) % Get a quick check which of models are loaded
clear names 
trim_data = [];
datvecx   = [];
datvecy   = [];
vecmag    = [];
for i = 1:length(tridx)
    variable     = 'depth averaged velocity';                       
    trim_data    = qpread(trim(tridx(i)),1,variable,'griddata',0,0,0);
    datvecx(:,:,i)      = (squeeze(mean(trim_data.XComp,1)))'; 
    datvecy(:,:,i)      = (squeeze(mean(trim_data.YComp,1)))';
    vecmag(:,:,i)  = (datvecx(:,:,i).^2+datvecy(:,:,i).^2).^(0.5);    
end
% Calculate differences between both model runs
dlt_datvecx = squeeze(datvecx(:,:,end)-datvecx(:,:,1)); 
dlt_datvecy = squeeze(datvecy(:,:,end)-datvecy(:,:,1));
dlt_vecmag  = squeeze(vecmag(:,:,end)-vecmag(:,:,1)); 
% Read Bottom depth: btdepth
trim_data = qpread(trim(tridx(end)),1,'bed level in water level points','griddata');
btdepth   = -(trim_data.Val(:,:))';

plotfig = figure('position',[50 50 1000 750]);
set(plotfig,'color','w');
% Plot Delta Flow Velocity Map Delta
sub(1) = subplot(3,2,[1,3,5]);
    % Plot Delta Magnitude
pcolor(X,Y,dlt_vecmag)
shading flat
shading interp
cb = colorbar;
set(cb,'Location','southoutside')
caxis(colax)
hold on
    % Plot Bottom Depth
contour(X,Y,btdepth,[-5:1.5:5],'ShowText','off','Color',[0.3 0.3 0.3])
    % Plot Diversion Area
plot([crdfig.x1div crdfig.x1div],[crdfig.y1div crdfig.y2div],'r')
plot([crdfig.x2div crdfig.x2div],[crdfig.y1div crdfig.y2div],'r')
plot([crdfig.x1div crdfig.x2div],[crdfig.y1div crdfig.y1div],'r')
plot([crdfig.x1div crdfig.x2div],[crdfig.y2div crdfig.y2div],'r')
    % Plot Polder Area
plot(Xr,Yr,'y','LineWidth',2)
    % Plot Flow velocities
quiver_thick(X(1:2:end,1:2:end),...
    Y(1:2:end,1:2:end),...
    dlt_datvecx(1:2:end,1:2:end),...
    dlt_datvecy(1:2:end,1:2:end),...
    'arrowhead_length',0.3,...
    'plot_color','y',...
    'uniform_length',1,...
    'colormap',bluewhitered(256),...
    'scaling',scale);
xlabel('Easting (km)')
ylabel('Northing (km)')
ylabel(cb,'Velocity difference (m/s)','fontsize',12)
axis equal
xlim([crdfig.x1 crdfig.x2]);
ylim([crdfig.y1 crdfig.y2]);
title_h = title(['\Delta Mean Flow Velocities (' nmmodel{tridx(end)} ' - ' nmmodel{tridx(1)} ')']);
title_x = crdfig.x2-0.15*(crdfig.x2-crdfig.x1);
title_y = crdfig.y2+0.05*(crdfig.y2-crdfig.y1);
set(title_h,'Position',[title_x title_y 0],'fontsize',12);
% Plot Flow Velocity Map Diversion
sub(2) = subplot(3,2,[4 6]);
    % Plot Delta Magnitude
pcolor(X,Y,dlt_vecmag)
shading flat
shading interp
caxis(colax)
hold on
    % Plot Bottom Depth
contour(X,Y,btdepth,[-5:1.5:5],'ShowText','on','Color',[0.3 0.3 0.3])
    % Plot Polder Area
plot(Xr,Yr,'y','LineWidth',2)
    % Plot Flow velocities
dlt_datvecx2 = dlt_datvecx;
dlt_datvecy2 = dlt_datvecy;
dlt_datvecx2(~in) = NaN;
dlt_datvecy2(~in) = NaN;
quiver_thick(X(1:2:end,1:2:end),...
    Y(1:2:end,1:2:end),...
    dlt_datvecx2(1:2:end,1:2:end),...
    dlt_datvecy2(1:2:end,1:2:end),...
    'arrowhead_length',0.3,...
    'plot_color','y',...
    'uniform_length',1,...
    'colormap',bluewhitered(256),...
    'scaling',scale);
    % Plot Scale Vector
x_vec = crdfig.x1div+0.025*(crdfig.x2div-crdfig.x1div);
x_txt = crdfig.x1div+0.025*(crdfig.x2div-crdfig.x1div);
y_vec = crdfig.y1div+0.05*(crdfig.y2div-crdfig.y1div);
y_txt = crdfig.y1div+0.025*(crdfig.y2div-crdfig.y1div);
quiver_thick(x_vec,y_vec,...
    u_vec,v_vec,...
    'arrowhead_length',0.3,...
    'plot_color','y',...
    'uniform_length',1,...
    'colormap',bluewhitered(256),...
    'scaling',scale);
text(x_txt,y_txt,[num2str(u_vec) ' m/s'],'color','k');
    % Plot Cross-sections
plot(coordcrs(:,1),coordcrs(:,2),'r')
set(gca,'xticklabel',[],'yticklabel',[])
axis equal
xlim([crdfig.x1div crdfig.x2div]);
ylim([crdfig.y1div crdfig.y2div])
set(sub(1),'position',[0.05 0.20 0.55 0.68]);  % Control position plot velocities Delta
set(sub(2),'position',[0.57 0.20 0.45 0.50]);  % Control position plot velocities Diversion

cd([fldr.figures fldr.nmdir 'Delta_VelocityMaps\']);
print(['DeltaMeanVelocityMap_' nmmodel{tridx(end)} '-' nmmodel{tridx(1)}],'-dpng');
cd(fldr.matlab);

%% 2 Delta Mean Sediment Transport and Sediment Accumulation

% Select model runs (never select more than 2 runs)
div_s = {'C'}; % Select from 'REF', 'A', 'B', or 'C'
q_s   = [290]; % Select only 1 discharge
sl_s  = [0 1]; % Select from '0' or '1'
% Select sediment class
sed_s = 'Sum Fractions' % Choose from 'Sum Fractions', 'Fine Sand', 'Medium Silt', or 'Fine Silt'
colax = [-0.05 0.05];

% Adjust size scale vector and color-axis based on discharge class 
switch q_s
    case {43}
        % Scale Vector Sediment Transport
        u_vec = 1E-6;   % Magnitude (m^3/s/m)
        v_vec = 0;
        scale = 1E5;    % Scaling
        % Multiplication factor for annual sedimentation
        m_factor = (0.5*365.25)/14.5  % multiplication factor = Days per year/Days per spring-neap cycle = Spring-neap cycles per year
    case {110}
        % Scale Vector Sediment Transport
        u_vec = 1E-5; % Magnitude (m^3/s/m)
        v_vec = 0;
        scale = 1E4;   % Scaling
        % Multiplication factor for annual sedimentation
        m_factor = (0.05*365.25)/14.5  % multiplication factor = Days per year/Days per sing-neap cycle        
    case {290}
        % Scale Vector Sediment Transport
        u_vec = 1E-4; % Magnitude (m^3/s/m)
        v_vec = 0;
        scale = 1E3
        % Multiplication factor for annual sedimentation
        m_factor = 0.5/14.5  % multiplication factor = Days per year/Days per spring-neap cycle 
    case {462}
        u_vec = 1E-3; % Magnitude (m^3/s/m)
        v_vec = 0;
        scale = 1E2;   % Scaling
        % Multiplication factor for annual sedimentation
        m_factor = 0.2/14.5  % multiplication factor = Days per year/Days per spring-neap cycle
end
% Select model runs
tridx = [];
for i = 1:length(div_s) 
    for j = 1:length(q_s)
        for k = 1:length(sl_s)
        tridx(i,j,k)  = find(strcmpi([div_s{i} '' sprintf('%d',find(q_s(j) == q)) '-sl' sprintf('%d',sl_s(k))],nmmodel));   
        end
    end
end
tridx = sort(tridx(:));
names = nmmodel(tridx) % Get a quick check which of models are loaded
clear names 
% Domains
    % Delta
crdfig.x1 = 5.20E2;
crdfig.x2 = 5.24E2;
crdfig.y1 = 52.13E2;
crdfig.y2 = 52.17E2;
    % Diversion
crdfig.x1div = 5.216E2;
crdfig.x2div = 5.225E2;
crdfig.y1div = 52.138E2;
crdfig.y2div = 52.147E2;
% Read sediment transport
    % Pre-allocate matrix size
S_ex    = zeros(length(tridx),size(X,1),size(X,2));
S_sx    = zeros(length(tridx),size(X,1),size(X,2));
S_ey    = zeros(length(tridx),size(X,1),size(X,2));
S_sy    = zeros(length(tridx),size(X,1),size(X,2));
datvecx = zeros(length(tridx),size(X,1),size(X,2));
datvecy = zeros(length(tridx),size(X,1),size(X,2));
sedacc  = zeros(length(tridx),(size(X,2)-1),(size(X,1)-1)); % Available mass of sediment is one cell smaller in X and Y
trim_data = [];
    % Read Transport and Accumulation Data 
for i = 1:length(tridx)
    cd([fldr.d3d fldr.nmdir fldr.model{tridx(i)}])
    switch sed_s
        case {'Sum Fractions'}
            trim_data    = load('meantotaltransport_sumfractions.mat');
            trim_data2   = load('availablemass_sumfractions.mat');
            rho          = 800; % Sediment Density (kg/m3)
        case {'Fine Sand'}
            trim_data    = load('meantotaltransport_finesand.mat');
            trim_data2   = load('availablemass_finesand.mat');
            rho          = 1200; % Sediment Density (kg/m3)
        case {'Medium Silt'}
            trim_data    = load('meantotaltransport_mediumsilt.mat');
            trim_data2   = load('availablemass_mediumsilt.mat');
            rho          = 600; % Sediment Density (kg/m3)
        case {'Fine Silt'}
            trim_data    = load('meantotaltransport_finesilt.mat');
            trim_data2   = load('availablemass_finesilt.mat');
            rho          = 400; % Sediment Density (kg/m3)
    end
trim_data.data.XComp = permute(trim_data.data.XComp,[3 2 1]); % Mean Transport in X (m3/s/m)
trim_data.data.YComp = permute(trim_data.data.YComp,[3 2 1]); % Mean Transport in Y (m3/s/m)
    % Calculate Mean Transport by correcting for spin-up
dt    = 60*60*24; % Conversion factor (days/sec) or (d/s)
t_run = 14.5;     % Simulation time without spin-up (d)
t_tot = t_run + spinup;
S_ex(i,:,:)    = squeeze((trim_data.data.XComp(:,:,end))).*dt*t_tot; % Total transport (m3/m) - run with spin-up
S_sx(i,:,:)    = squeeze((trim_data.data.XComp(:,:,1))).*dt*spinup;  % Total transport (m3/m) - only spin-up
datvecx(i,:,:) = (S_ex(i,:,:) - S_sx(i,:,:))/(t_run*dt);                         % Mean Transport in X (m3/s/m) - run corrected for spin-up
S_ey(i,:,:)    = squeeze((trim_data.data.YComp(:,:,end))).*dt*t_tot; % Total transport (m3/m) - run with spin-up
S_sy(i,:,:)    = squeeze((trim_data.data.YComp(:,:,1))).*dt*spinup;  % Total transport (m3/m) - only spin-up
datvecy(i,:,:) = (S_ey(i,:,:) - S_sy(i,:,:))/(t_run*dt);  
clear S_ex S_sx S_ey S_sy
    % Calculate Sediment accumulation
trim_data2.data.X = squeeze(trim_data2.data.X(:,2:end,2:end,1));
trim_data2.data.Y = squeeze(trim_data2.data.Y(:,2:end,2:end,1));
Xs             = (trim_data2.data.X./1000)'; % Convert from m to km
Ys             = (trim_data2.data.Y./1000)'; % Convert from m to km
sedacc(i,:,:)  = squeeze((trim_data2.data.Val(end,:,:)-trim_data2.data.Val(1,:,:)))./rho; % Net sediment accumulation in (m)
end
cd([fldr.matlab])
sedacc  = permute(sedacc,[1 3 2]);
    % Calculate differences between both model runs (Delta)
dlt_datvecx = squeeze(datvecx(end,:,:)-datvecx(1,:,:));
dlt_datvecy = squeeze(datvecy(end,:,:)-datvecy(1,:,:));
dlt_sedacc  = squeeze(sedacc(end,:,:)-sedacc(1,:,:));
% Diffrence in Potential Annual Accumulation (Correct for occurence of a discharge class)
dlt_sedacc  = m_factor.*dlt_sedacc; 
% Read Bottom depth: btdepth
trim_data = qpread(trim(tridx(end)),1,'bed level in water level points','griddata');
btdepth   = -(trim_data.Val(:,:))';

plotfig = figure('position',[50 50 1000 750]);
set(plotfig,'color','w');
% Plot Transport Map Delta
sub(1) = subplot(3,2,[1,3,5]);
    % Plot Sediment Accumulation
pcolor(Xs,Ys,dlt_sedacc);
shading flat
shading interp
cb = colorbar;
set(cb,'Location','southoutside')
caxis(colax)
hold on
    % Plot Bottom Depth
contour(X,Y,btdepth,[-5:1.5:5],'ShowText','off','Color',[0.3 0.3 0.3])
    % Plot Diversion Area
plot([crdfig.x1div crdfig.x1div],[crdfig.y1div crdfig.y2div],'r')
plot([crdfig.x2div crdfig.x2div],[crdfig.y1div crdfig.y2div],'r')
plot([crdfig.x1div crdfig.x2div],[crdfig.y1div crdfig.y1div],'r')
plot([crdfig.x1div crdfig.x2div],[crdfig.y2div crdfig.y2div],'r')
    % Plot Polder Area
plot(Xr,Yr,'y','LineWidth',1.5)
    % Plot Sediment Transport
quiver_thick(X(1:2:end,1:2:end),...
    Y(1:2:end,1:2:end),...
    dlt_datvecx(1:2:end,1:2:end),...
    dlt_datvecy(1:2:end,1:2:end),...
    'arrowhead_length',0.3,...
    'plot_color','y',...
    'uniform_length',1,...
    'colormap',bluewhitered(256),...
    'scaling',scale);
xlabel('Easting (km)')
ylabel('Northing (km)')
ylabel(cb,'Accumulation difference (m)','fontsize',12)
axis equal
xlim([crdfig.x1 crdfig.x2]);
ylim([crdfig.y1 crdfig.y2]);
title_h = title(['\Delta Mean Transport and Potential Annual Accumulation of ' sed_s ' (' nmmodel{tridx(end)} ' - ' nmmodel{tridx(1)} ')']);
title_x = crdfig.x2-0.15*(crdfig.x2-crdfig.x1);
title_y = crdfig.y2+0.05*(crdfig.y2-crdfig.y1);
set(title_h,'Position',[title_x title_y 0],'fontsize',12);
% Plot Transport Map Diversion
sub(2) = subplot(3,2,[4 6]);
    % Plot Sediment Accumulation
pcolor(Xs,Ys,dlt_sedacc);
shading flat
shading interp
caxis(colax)
hold on
    % Plot Bottom Depth
contour(X,Y,btdepth,[-5:1.5:5],'ShowText','on','Color',[0.3 0.3 0.3])
    % Plot Polder Area
plot(Xr,Yr,'y','LineWidth',1.5)
    % Plot Sediment Transport
dlt_datvecx2 = dlt_datvecx;
dlt_datvecy2 = dlt_datvecy;
dlt_datvecx2(~in) = NaN;
dlt_datvecy2(~in) = NaN;
quiver_thick(X(1:2:end,1:2:end),...
    Y(1:2:end,1:2:end),...
    dlt_datvecx2(1:2:end,1:2:end),...
    dlt_datvecy2(1:2:end,1:2:end),...
    'arrowhead_length',0.3,...
    'plot_color','y',...
    'uniform_length',1,...
    'colormap',bluewhitered(256),...
    'scaling',scale);
x_vec = crdfig.x1div+0.025*(crdfig.x2div-crdfig.x1div);
x_txt = crdfig.x1div+0.025*(crdfig.x2div-crdfig.x1div);
y_vec = crdfig.y1div+0.05*(crdfig.y2div-crdfig.y1div);
y_txt = crdfig.y1div+0.025*(crdfig.y2div-crdfig.y1div);
    % Plot Scale Vector
quiver_thick(x_vec,y_vec,...
    u_vec,v_vec,...
    'arrowhead_length',0.3,...
    'plot_color','y',...
    'uniform_length',1,...
    'colormap',bluewhitered(256),...
    'scaling',scale);
text(x_txt,y_txt,[num2str(u_vec) ' m3/s/m'],'color','k');
set(gca,'xticklabel',[],'yticklabel',[])
axis equal
xlim([crdfig.x1div crdfig.x2div]);
ylim([crdfig.y1div crdfig.y2div])
set(sub(1),'position',[0.05 0.20 0.55 0.68]);   % Control position plot velocities Delta
set(sub(2),'position',[0.57 0.20 0.45 0.5]);   % Control position plot velocities Diversion

cd([fldr.figures fldr.nmdir 'Delta_TransportMaps\']);
print(['DeltaMeanTransportAccumulationMap_' nmmodel{tridx(end)} '-' nmmodel{tridx(1)}],'-dpng');
cd(fldr.matlab);

%% 3 Restoration Impact due to SLR: Delta Mean Sediment Transport and Sediment Accumulation (hard-coded)

% Select model runs (never select more than 2 runs)
q_s   = [110]; % Select only 1 discharge
div_s = {'REF','C'}; % Select from 'REF', 'A', 'B', or 'C'
sl_s  = [0 1];       
% Select sediment class
sed_s = 'Sum Fractions' % Choose from 'Sum Fractions', 'Fine Sand', 'Medium Silt', or 'Fine Silt'
colax = [-0.05 0.05];

% Adjust size scale vector and color-axis based on discharge class 
switch q_s
    case {43}
        % Scale Vector Sediment Transport
        u_vec = 1E-6;   % Magnitude (m^3/s/m)
        v_vec = 0;
        scale = 5E4;    % Scaling
        % Multiplication factor for annual sedimentation
        m_factor = (0.5*365.25)/14.5  % multiplication factor = Days per year/Days per spring-neap cycle = Spring-neap cycles per year
    case {110}
        % Scale Vector Sediment Transport
        u_vec = 1E-6; % Magnitude (m^3/s/m)
        v_vec = 0;
        scale = 1E5;   % Scaling
        % Multiplication factor for annual sedimentation
        m_factor = (0.05*365.25)/14.5  % multiplication factor = Days per year/Days per sing-neap cycle        
    case {290}
        % Scale Vector Sediment Transport
        u_vec = 1E-5; % Magnitude (m^3/s/m)
        v_vec = 0;
        scale = 1E4;   % Scaling
        % Multiplication factor for annual sedimentation
        m_factor = 0.5/14.5  % multiplication factor = Days per year/Days per spring-neap cycle 
    case {462}
        u_vec = 1E-4; % Magnitude (m^3/s/m)
        v_vec = 0;
        scale = 1E3;   % Scaling
        % Multiplication factor for annual sedimentation
        m_factor = 0.2/14.5  % multiplication factor = Days per year/Days per spring-neap cycle
end
% Select model runs
tridx = [];
for i = 1:length(div_s) 
    for j = 1:length(q_s)
        for k = 1:length(sl_s)
        tridx(i,j,k)  = find(strcmpi([div_s{i} '' sprintf('%d',find(q_s(j) == q)) '-sl' sprintf('%d',sl_s(k))],nmmodel));   
        end
    end
end
tridx = sort(tridx(:));
names = nmmodel(tridx) % Get a quick check which of models are loaded
clear names 
% Domains
    % Delta
crdfig.x1 = 5.20E2;
crdfig.x2 = 5.24E2;
crdfig.y1 = 52.13E2;
crdfig.y2 = 52.17E2;
    % Diversion
crdfig.x1div = 5.216E2;
crdfig.x2div = 5.225E2;
crdfig.y1div = 52.138E2;
crdfig.y2div = 52.147E2;
% Read sediment transport
    % Pre-allocate matrix size
S_ex    = zeros(length(tridx),size(X,1),size(X,2));
S_sx    = zeros(length(tridx),size(X,1),size(X,2));
S_ey    = zeros(length(tridx),size(X,1),size(X,2));
S_sy    = zeros(length(tridx),size(X,1),size(X,2));
datvecx = zeros(length(tridx),size(X,1),size(X,2));
datvecy = zeros(length(tridx),size(X,1),size(X,2));
sedacc  = zeros(length(tridx),(size(X,2)-1),(size(X,1)-1)); % Available mass of sediment is one cell smaller in X and Y
trim_data = [];
    % Read Transport and Accumulation Data 
for i = 1:length(tridx)
    cd([fldr.d3d fldr.nmdir fldr.model{tridx(i)}])
    switch sed_s
        case {'Sum Fractions'}
            trim_data    = load('meantotaltransport_sumfractions.mat');
            trim_data2   = load('availablemass_sumfractions.mat');
            rho          = 800; % Sediment Density (kg/m3)
        case {'Fine Sand'}
            trim_data    = load('meantotaltransport_finesand.mat');
            trim_data2   = load('availablemass_finesand.mat');
            rho          = 1200; % Sediment Density (kg/m3)
        case {'Medium Silt'}
            trim_data    = load('meantotaltransport_mediumsilt.mat');
            trim_data2   = load('availablemass_mediumsilt.mat');
            rho          = 600; % Sediment Density (kg/m3)
        case {'Fine Silt'}
            trim_data    = load('meantotaltransport_finesilt.mat');
            trim_data2   = load('availablemass_finesilt.mat');
            rho          = 400; % Sediment Density (kg/m3)
    end
trim_data.data.XComp = permute(trim_data.data.XComp,[3 2 1]); % Mean Transport in X (m3/s/m)
trim_data.data.YComp = permute(trim_data.data.YComp,[3 2 1]); % Mean Transport in Y (m3/s/m)
    % Calculate Mean Transport by correcting for spin-up
dt    = 60*60*24; % Conversion factor (days/sec) or (d/s)
t_run = 14.5;     % Simulation time without spin-up (d)
t_tot = t_run + spinup;
S_ex(i,:,:)    = squeeze((trim_data.data.XComp(:,:,end))).*dt*t_tot; % Total transport (m3/m) - run with spin-up
S_sx(i,:,:)    = squeeze((trim_data.data.XComp(:,:,1))).*dt*spinup;  % Total transport (m3/m) - only spin-up
datvecx(i,:,:) = (S_ex(i,:,:) - S_sx(i,:,:))/(t_run*dt);                         % Mean Transport in X (m3/s/m) - run corrected for spin-up
S_ey(i,:,:)    = squeeze((trim_data.data.YComp(:,:,end))).*dt*t_tot; % Total transport (m3/m) - run with spin-up
S_sy(i,:,:)    = squeeze((trim_data.data.YComp(:,:,1))).*dt*spinup;  % Total transport (m3/m) - only spin-up
datvecy(i,:,:) = (S_ey(i,:,:) - S_sy(i,:,:))/(t_run*dt);  
clear S_ex S_sx S_ey S_sy
    % Calculate Sediment accumulation
trim_data2.data.X = squeeze(trim_data2.data.X(:,2:end,2:end,1));
trim_data2.data.Y = squeeze(trim_data2.data.Y(:,2:end,2:end,1));
Xs             = (trim_data2.data.X./1000)'; % Convert from m to km
Ys             = (trim_data2.data.Y./1000)'; % Convert from m to km
sedacc(i,:,:)  = squeeze((trim_data2.data.Val(end,:,:)-trim_data2.data.Val(1,:,:)))./rho; % Net sediment accumulation in (m)
end
cd([fldr.matlab])
sedacc  = permute(sedacc,[1 3 2]);
    % Calculate differences between both model runs (Delta)
dlt_datvecx_sl0 = squeeze(datvecx(2,:,:)-datvecx(1,:,:));
dlt_datvecy_sl0 = squeeze(datvecy(2,:,:)-datvecy(1,:,:));
dlt_sedacc_sl0  = squeeze(sedacc(2,:,:)-sedacc(1,:,:));
dlt_sedacc_sl0  = m_factor.*dlt_sedacc_sl0; 
    % Calculate differences between both model runs (Delta)
dlt_datvecx_sl1 = squeeze(datvecx(4,:,:)-datvecx(3,:,:));
dlt_datvecy_sl1 = squeeze(datvecy(4,:,:)-datvecy(3,:,:));
dlt_sedacc_sl1  = squeeze(sedacc(4,:,:)-sedacc(3,:,:));
dlt_sedacc_sl1  = m_factor.*dlt_sedacc_sl1;
% Calculate the difference of the difference, i.e. the combined impact of
% restoration (C-REF) and SLR (SL1-SL0)
dlt_dlt_datvecx = dlt_datvecx_sl1 - dlt_datvecx_sl0;
dlt_dlt_datvecy = dlt_datvecy_sl1 - dlt_datvecy_sl0;
dlt_dlt_sedacc  = dlt_sedacc_sl1 - dlt_sedacc_sl0;
% Read Bottom depth: btdepth
trim_data = qpread(trim(tridx(end)),1,'bed level in water level points','griddata');
btdepth   = -(trim_data.Val(:,:))';

plotfig = figure('position',[50 50 1000 750]);
set(plotfig,'color','w');
% Plot Transport Map Delta
sub(1) = subplot(3,2,[1,3,5]);
    % Plot Sediment Accumulation
pcolor(Xs,Ys,dlt_dlt_sedacc);
shading flat
shading interp
cb = colorbar;
set(cb,'Location','southoutside')
caxis(colax)
hold on
    % Plot Bottom Depth
contour(X,Y,btdepth,[-5:1.5:5],'ShowText','off','Color',[0.3 0.3 0.3])
    % Plot Diversion Area
plot([crdfig.x1div crdfig.x1div],[crdfig.y1div crdfig.y2div],'r')
plot([crdfig.x2div crdfig.x2div],[crdfig.y1div crdfig.y2div],'r')
plot([crdfig.x1div crdfig.x2div],[crdfig.y1div crdfig.y1div],'r')
plot([crdfig.x1div crdfig.x2div],[crdfig.y2div crdfig.y2div],'r')
    % Plot Polder Area
plot(Xr,Yr,'y','LineWidth',1.5)
    % Plot Sediment Transport
quiver_thick(X(1:2:end,1:2:end),...
    Y(1:2:end,1:2:end),...
    dlt_dlt_datvecx(1:2:end,1:2:end),...
    dlt_dlt_datvecy(1:2:end,1:2:end),...
    'arrowhead_length',0.3,...
    'plot_color','y',...
    'uniform_length',1,...
    'colormap',bluewhitered(256),...
    'scaling',scale);
xlabel('Easting (km)')
ylabel('Northing (km)')
ylabel(cb,'Accumulation difference (m)','fontsize',12)
axis equal
xlim([crdfig.x1 crdfig.x2]);
ylim([crdfig.y1 crdfig.y2]);
title_h = title(['\Delta Mean Transport and Potential Annual Accumulation of ' sed_s ' {(' nmmodel{tridx(4)} ' - ' nmmodel{tridx(3)} ') - (' nmmodel{tridx(2)} ' - ' nmmodel{tridx(1)} ')}']);
title_x = crdfig.x2-0.15*(crdfig.x2-crdfig.x1);
title_y = crdfig.y2+0.05*(crdfig.y2-crdfig.y1);
set(title_h,'Position',[title_x title_y 0],'fontsize',12);
% Plot Transport Map Diversion
sub(2) = subplot(3,2,[4 6]);
    % Plot Sediment Accumulation
pcolor(Xs,Ys,dlt_dlt_sedacc);
shading flat
shading interp
caxis(colax)
hold on
    % Plot Bottom Depth
contour(X,Y,btdepth,[-5:1.5:5],'ShowText','on','Color',[0.3 0.3 0.3])
    % Plot Polder Area
plot(Xr,Yr,'y','LineWidth',1.5)
    % Plot Sediment Transport
dlt_dlt_datvecx2 = dlt_dlt_datvecx;
dlt_dlt_datvecy2 = dlt_dlt_datvecy;
dlt_dlt_datvecx2(~in) = NaN;
dlt_dlt_datvecy2(~in) = NaN;
quiver_thick(X(1:2:end,1:2:end),...
    Y(1:2:end,1:2:end),...
    dlt_dlt_datvecx2(1:2:end,1:2:end),...
    dlt_dlt_datvecy2(1:2:end,1:2:end),...
    'arrowhead_length',0.3,...
    'plot_color','y',...
    'uniform_length',1,...
    'colormap',bluewhitered(256),...
    'scaling',scale);
x_vec = crdfig.x1div+0.025*(crdfig.x2div-crdfig.x1div);
x_txt = crdfig.x1div+0.025*(crdfig.x2div-crdfig.x1div);
y_vec = crdfig.y1div+0.05*(crdfig.y2div-crdfig.y1div);
y_txt = crdfig.y1div+0.025*(crdfig.y2div-crdfig.y1div);
    % Plot Scale Vector
quiver_thick(x_vec,y_vec,...
    u_vec,v_vec,...
    'arrowhead_length',0.3,...
    'plot_color','y',...
    'uniform_length',1,...
    'colormap',bluewhitered(256),...
    'scaling',scale);
text(x_txt,y_txt,[num2str(u_vec) ' m3/s/m'],'color','k');
set(gca,'xticklabel',[],'yticklabel',[])
axis equal
xlim([crdfig.x1div crdfig.x2div]);
ylim([crdfig.y1div crdfig.y2div])
set(sub(1),'position',[0.05 0.20 0.55 0.68]);   % Control position plot velocities Delta
set(sub(2),'position',[0.57 0.20 0.45 0.5]);   % Control position plot velocities Diversion

cd([fldr.figures fldr.nmdir 'Delta_TransportMaps\']);
print(['Delta_DeltaMeanTransportAccumulationMap_SLR_' div_s{end} '-' div_s{1} '_'  num2str(q_s)],'-dpng');
cd(fldr.matlab);