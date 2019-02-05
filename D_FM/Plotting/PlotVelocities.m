clearvars
% Load in the grid and subsample to be just bbay
ldb=landboundary('read',['C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\FM\LandBoundaryFiles\',...
    'WA_coastline.ldb']);
dname='E:\DelftFM\AndrewRuns\RoughnessTesting\Manning\r_316\';
fname    = 'ps_2d_map.nc';
grd         = dflowfm.readNet([dname,fname]);
% ---------------- THIS IS ALL FOR GETTING TRIANGLES TO PLOT---------------
xg = grd.node.x(grd.edge.NetLink);
yg = grd.node.y(grd.edge.NetLink);
xg(3,:)=NaN;
yg(3,:)=NaN;
inds = yg(1,:) <= 48.83 & xg(1,:) >= -122.77 & yg(1,:) >= 48.6 & xg(1,:) <= -122.35;
xg(:,~inds) = [];
yg(:,~inds) = [];
plotting = 0;
if plotting
    f=figure;
    set(f,'renderer','zbuffer','units','inches',...
        'position',[0 0 9 7],...
        'paperpositionmode','au')
    axesm('mercator')
    lh=plotm(yg(:),xg(:),'-');
    hold on
    plotm(ldb(:,2),ldb(:,1),'k-')
end
tim = nc_cf_time([dname,fname]);
info = ncinfo([dname fname]);
% -------------------------------------------------------------------------
inds = grd.face.FlowElem_y <= 48.83 & grd.face.FlowElem_x >= -122.77 & grd.face.FlowElem_y >= 48.6 & grd.face.FlowElem_x <= -122.35;
yg = grd.face.FlowElem_y(inds);
xg = grd.face.FlowElem_x(inds);

% Load in cherry point water level
load('E:\DelftFM\data\wl\coops_2017\navd88\NOAA_9449424_Cherry_Point_navd88.mat\');
tvec = datenum(2017,1,1):1/24:datenum(2017,12,31);
wld.wl = interp1(wld.time,wld.WL_VALUE,tvec); % Convert to hourly
tin1 = find(tvec == tim(1)); % Find where time overlaps with model time
tin2 = find(tvec == tim(end));
wld.time = tvec(tin1:tin2);
wld.wl = wld.wl(tin1:tin2);



%% Load in velocities
% Note: Row and Column is in first argument of ncread and then how many to go in that direction is the second.
set(0,'defaultAxesFontSize',8)
% Dimensions for plot
p_left = .05;
p_right = .05;
p_top = .02;
p_bot = .08;
p_spacing = 0.04;
p_big = .05;
p_wid = (1-p_right-p_left);
p_height = (1-p_top-p_bot)/5;

% Time frame of concern
tstart = datenum(2017,2,7);
tend = datenum(2017,2,23);
is = find(wld.time == tstart);
ie = find(wld.time == tend);
tvec = wld.time(is):1/24:wld.time(ie);
vec = is:1:ie;

v = VideoWriter('BellinghamBay_02072017_02232017','MPEG-4');
v.FrameRate = 3;
v.Quality = 75;
open(v)

% Load in Basemap
IM = load('C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\GeoTiffs\Bbay\spherical\bbay_chips.mat');

for ii = 1:length(vec)
    clf
    vx = ncread([dname,fname],'ucx',[1 vec(ii)],[394013 1]); % Velocities
    vy = ncread([dname,fname],'ucy',[1 vec(ii)],[394013 1]);
    vx = vx(inds);
    vy = vy(inds);
    axes('position',[p_left p_bot+p_height p_wid p_height*4]);
    imagesc(IM.xm,IM.ym,IM.im)
    set(gca,'ydir','normal')
    hold on
    xlabel('Degrees Longitude')
    ylabel('Degrees Latitude')
    q = quiver(xg(1:2:end),yg(1:2:end),vx(1:2:end)',vy(1:2:end)','Color','w','AutoScaleFactor',2.5);
    
    % Add tide level on bottom with moving line to show time
    axes('position',[p_left p_bot-p_spacing p_wid p_height]);
    plot(wld.time(is:ie),wld.wl(is:ie))
    hold on
    line([wld.time(vec(ii)) wld.time(vec(ii))],[-1 4],'Color','r')
    datetick()
    ylabel('Tide Level [m]')
    xlabel('Time')
    
    %     pause(.8)
    frame = getframe(gcf);
    writeVideo(v,frame);
    
end
close(v);

%% Plot of averaged flow over spring/neap cycle
clf
% Load in Basemap
IM = load('C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\GeoTiffs\Bbay\spherical\bbay_chips_r2.mat');
I = imagesc(IM.xm,IM.ym,IM.im);
set(gca,'ydir','normal')
hold on
% BBay Limits
% xlim([-122.653 -122.452])
ylim([48.69 max(IM.ym)])

% avgX = zeros(size(vx));
% avgY = zeros(size(vy));
% for ii = 1:length(vec)
%     vx = ncread([dname,fname],'ucx',[1 vec(ii)],[394013 1]); % Velocities
%     vy = ncread([dname,fname],'ucy',[1 vec(ii)],[394013 1]);
%     avgX = avgX + vx(inds);
%     avgY = avgY + vy(inds);
% end
    
% q = quiver(xg(1:2:end),yg(1:2:end),avgX(1:2:end)',avgY(1:2:end)','Color','w','AutoScaleFactor',5);
q = quiver(xg,yg,avgX',avgY','Color','w','AutoScale','off');

xlabel('Degrees Longitude')
ylabel('Degrees Latitude')
set(gca,'FontSize',12)
printFig(gcf,'Bay_AVGFlow_02072017_02232017',[15 12],'png',300)
%% Stefans Code - Miss you stefan 









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


