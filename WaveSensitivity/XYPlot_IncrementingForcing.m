clearvars
addpath C:\Functions_Matlab
addpath C:\Functions_Matlab\mapping\kml

kml = dir('C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\WaveAnalysis\Sensitivity\KML_Points\*.kml');
kmlfol = 'C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\WaveAnalysis\Sensitivity\KML_Points\';
format = [1 4 2 3];
for ii = 1:length(kml)
    t = kml2struct([kmlfol kml(format(ii)).name]);
    K(ii).lon = t.Lon;
    K(ii).lat = t.Lat;
    K(ii).name = t.Name;
    [K(ii).x,K(ii).y] = deg2utm(t.Lat,t.Lon);
end



return
%% Sensitivity to speed 
direction = 'south';
clf
switch direction   % 'Coming from the....' notion 
    case 'north'
        dirWant = 0;
        dir2use = wrap2360(90-(dirWant - 180));
    case 'south'
        dirWant = 180;
        dir2use = wrap2360(90-(dirWant - 180));
        fname1 = sprintf('southern_t00_s20_d%d.mat',dir2use);
        fname2 = 'SpatialFromSouth.mat';
    otherwise
        disp('Error with direction')
end
b = parula(5);
spd = 5;
transform = 0;
for ii = 1:length(K)
    hsig = zeros(1,10);
    for jj = 1:length(hsig)
        spd2load = round(spd);
        fname = sprintf('southern_t200_s%d_d%d.mat',spd2load,dir2use);
        D = load(['C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\Hindcast\Southern_LUT_Data\' fname]);
        [xg,yg] = findNearestGridPoint(D.X,D.Y,K(ii).lon,K(ii).lat); % Find closest wave grid point
        hsig(jj) = D.hs(xg,yg); % Grab wave point data
        if transform
            hsig(jj) = D.hs(xg,yg)^2;
        end
        spd = spd+2.5; %Increase Speed
    end
    x_axis = 5:2.5:5+(2.5*length(hsig)-1);
    ll(ii) = plot(x_axis,hsig,'o','MarkerSize',9,'MarkerFaceColor',b(ii,:),'MarkerEdgeColor',b(ii,:));
    hold on 
    line(x_axis,hsig,'Color',b(ii,:));
    hold on
    clear hsig 
    spd = 5;
end
grid on 
box on 
xlabel('Wind Speed [m/s]','FontSize',12)
ylabel('Hsig [m]','FontSize',12)
ax = gca;
ax.FontSize = 12;
lgd = legend(ll,arrayfun(@(x)(x.name),K,'un',0),...
    'interpreter','none',...
    'location','northwest');
lgd.FontSize = 12;
fname = sprintf('SpeedSensitivity_WavesAtPointsFrom%s',upper(direction));
if transform
    fname = sprintf('SpeedSensitivity_WavesAtPointsFrom%s_SquaredTransform',upper(direction));
end
printFig(gcf,fname,[12 11],'png',300)
%% Sensitivity to Direction
direction = 'south';
clf
switch direction   % 'Coming from the....' notion 
    case 'north'
        dirWant = 0;
        dir2use = wrap2360(90-(dirWant - 180));
    case 'south'
        dirWant = 180;
        dir2use = wrap2360(90-(dirWant - 180));
        fname1 = sprintf('southern_t00_s20_d%d.mat',dir2use);
        fname2 = 'SpatialFromSouth.mat';
    otherwise
        disp('Error with direction')
end
b = parula(5);
dirWant = 0;
for ii = 1:length(K)
    hsig = zeros(1,24);
    for jj = 1:length(hsig)
        spd2load = 15;
        dir2use = wrap2360(90-(dirWant - 180));
        fname = sprintf('southern_t200_s%d_d%d.mat',spd2load,dir2use);
        D = load(['C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\Hindcast\Southern_LUT_Data\' fname]);
        [xg,yg] = findNearestGridPoint(D.X,D.Y,K(ii).lon,K(ii).lat); % Find closest wave grid point
        hsig(jj) = D.hs(xg,yg); % Grab wave point data
        dirWant = dirWant + 15; %Increase Speed
    end
    x_axis = 0:15:345;
    ll(ii) = plot(x_axis,hsig,'o','MarkerSize',9,'MarkerFaceColor',b(ii,:),'MarkerEdgeColor',b(ii,:));
    hold on 
    line(x_axis,hsig,'Color',b(ii,:));
    hold on
    clear hsig 
    dirWant = 0;
end
grid on 
box on 
xlabel('Wind Direction [degrees]','FontSize',12)
ylabel('Hsig [m]','FontSize',12)
ax = gca;
ax.FontSize = 12;
lgd = legend(ll,arrayfun(@(x)(x.name),K,'un',0),...
    'interpreter','none',...
    'location','north');
lgd.FontSize = 12;
fname = sprintf('DirectionSensitivity_WavesAtPoints_Speed_%d',spd2load);
printFig(gcf,fname,[12 11],'png',300)
%% Sensitivity to Tide Level 

direction = 'south';
clf
switch direction   % 'Coming from the....' notion 
    case 'north'
        dirWant = 0;
        dir2use = wrap2360(90-(dirWant - 180));
    case 'south'
        dirWant = 180;
        dir2use = wrap2360(90-(dirWant - 180));
        fname1 = sprintf('southern_t00_s20_d%d.mat',dir2use);
        fname2 = 'SpatialFromSouth.mat';
    otherwise
        disp('Error with direction')
end
b = parula(5);
tide2load = -200;
for ii = 1:length(K)
    hsig = zeros(1,16);
    for jj = 1:length(hsig)
        spd2load = 15;    
        fname = sprintf('southern_t%i_s%d_d%d.mat',tide2load,spd2load,dir2use);
        if tide2load == 0
            tide2load = '00';
            fname = sprintf('southern_t%s_s%d_d%d.mat',tide2load,spd2load,dir2use);
            tide2load = 0;
        end
        D = load(['C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\Hindcast\Southern_LUT_Data\' fname]);
        [xg,yg] = findNearestGridPoint(D.X,D.Y,K(ii).lon,K(ii).lat); % Find closest wave grid point
        hsig(jj) = D.hs(xg,yg); % Grab wave point data
        tide2load = tide2load + 50; %Increase Speed
    end
    x_axis = -2:.5:5.5;
    ll(ii) = plot(x_axis,hsig,'o','MarkerSize',9,'MarkerFaceColor',b(ii,:),'MarkerEdgeColor',b(ii,:));
    hold on 
    line(x_axis,hsig,'Color',b(ii,:));
    hold on
    clear hsig 
    tide2load = -200;
end
grid on 
box on 
xlabel('Water Level [m]','FontSize',12)
ylabel('Hsig [m]','FontSize',12)
ax = gca;
ax.FontSize = 12;
lgd = legend(ll,arrayfun(@(x)(x.name),K,'un',0),...
    'interpreter','none',...
    'location','west');
lgd.FontSize = 12;
fname = sprintf('Waterlevel_Sensitivity_WavesAtPoints_Speed_%d_WindsFromThe_%s',spd2load,upper(direction));
printFig(gcf,fname,[12 11],'png',300)

%% Combined Figure of Direction and Wind Speed

% Specify Dimensions of Plot
p_left = .05;
p_right = .05;
p_top = .03;
p_bot = .07;
p_spacing = .01;
p_wid = ((1-p_right-p_left-p_spacing)/2);
p_height = (1-p_top-p_bot-p_spacing);

row = 0;
col = 0;
% -------------- Speed Sensitivity 
direction = 'north';
clf
switch direction   % 'Coming from the....' notion 
    case 'north'
        dirWant = 0;
        dir2use = wrap2360(90-(dirWant - 180));
    case 'south'
        dirWant = 180;
        dir2use = wrap2360(90-(dirWant - 180));
        fname1 = sprintf('southern_t00_s20_d%d.mat',dir2use);
        fname2 = 'SpatailFromSouth.mat';
    otherwise
        disp('Error with direction')
end
b = parula(5);
spd = 5;
axes('position',[p_left+col*(p_wid+p_spacing) p_bot+row*(p_height+p_spacing) p_wid p_height])
for ii = 1:length(K)
    hsig = zeros(1,10);
    for jj = 1:length(hsig)
        spd2load = round(spd);
        fname = sprintf('southern_t200_s%d_d%d.mat',spd2load,dir2use);
        D = load(['C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\Hindcast\Southern_LUT_Data\' fname]);
        [xg,yg] = findNearestGridPoint(D.X,D.Y,K(ii).lon,K(ii).lat); % Find closest wave grid point
        hsig(jj) = D.hs(xg,yg); % Grab wave point data
        spd = spd+2.5; %Increase Speed
    end
    x_axis = 5:2.5:5+(2.5*length(hsig)-1);
    ll(ii) = plot(x_axis,hsig,'o','MarkerSize',9,'MarkerFaceColor',b(ii,:),'MarkerEdgeColor',b(ii,:));
    hold on 
    line(x_axis,hsig,'Color',b(ii,:));
    hold on
    clear hsig 
    spd = 5;
end
grid on 
box on 
xlabel('Wind Speed [m/s]','FontSize',12)
ylabel('Hsig [m]','FontSize',12)
ax = gca;
ax.FontSize = 12;
lgd = legend(ll,arrayfun(@(x)(x.name),K,'un',0),...
    'interpreter','none',...
    'location','northwest');
lgd.FontSize = 12;

% -------------- Direction Sensitivity 
col = 1;
direction = 'south';
b = parula(5);
dirWant = 0;
axes('position',[p_left+col*(p_wid+p_spacing)+.028 p_bot+row*(p_height+p_spacing) p_wid p_height])
for ii = 1:length(K)
    hsig = zeros(1,24);
    for jj = 1:length(hsig)
        spd2load = 15;
        dir2use = wrap2360(90-(dirWant - 180));
        fname = sprintf('southern_t200_s%d_d%d.mat',spd2load,dir2use);
        D = load(['C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\Hindcast\Southern_LUT_Data\' fname]);
        [xg,yg] = findNearestGridPoint(D.X,D.Y,K(ii).lon,K(ii).lat); % Find closest wave grid point
        hsig(jj) = D.hs(xg,yg); % Grab wave point data
        dirWant = dirWant + 15; %Increase Speed
    end
    x_axis = 0:15:345;
    ll(ii) = plot(x_axis,hsig,'o','MarkerSize',9,'MarkerFaceColor',b(ii,:),'MarkerEdgeColor',b(ii,:));
    hold on 
    line(x_axis,hsig,'Color',b(ii,:));
    hold on
    clear hsig 
    dirWant = 0;
end
grid on 
box on 
xlabel('Wind Direction [degrees]','FontSize',12)
ylabel(' ','FontSize',12)
ax = gca;
ax.FontSize = 12;
% lgd = legend(ll,arrayfun(@(x)(x.name),K,'un',0),...
%     'interpreter','none',...
%     'location','north');
% lgd.FontSize = 12;

fname = sprintf('WaveSensitivty_Speed_Direction_SpeedOf%d',spd2load);
printFig(gcf,fname,[11 8.5],'png',300)

%% Combined Figure of Direction and Wind Speed and Water level sensitify of Locations 

% Specify Dimensions of Plot
p_left = .05;
p_right = .05;
p_top = .03;
p_bot = .07;
p_spacing = .01;
p_wid = ((1-p_right-p_left-p_spacing)/3.3);
p_height = (1-p_top-p_bot-p_spacing);

row = 0;
col = 0;
% -------------- Speed Sensitivity 
direction = 'north';
clf
switch direction   % 'Coming from the....' notion 
    case 'north'
        dirWant = 0;
        dir2use = wrap2360(90-(dirWant - 180));
    case 'south'
        dirWant = 180;
        dir2use = wrap2360(90-(dirWant - 180));
        fname1 = sprintf('southern_t00_s20_d%d.mat',dir2use);
        fname2 = 'SpatailFromSouth.mat';
    otherwise
        disp('Error with direction')
end
b = parula(5);
spd = 5;
axes('position',[p_left+col*(p_wid+p_spacing) p_bot+row*(p_height+p_spacing) p_wid p_height])
for ii = 1:length(K)
    hsig = zeros(1,10);
    for jj = 1:length(hsig)
        spd2load = round(spd);
        fname = sprintf('southern_t200_s%d_d%d.mat',spd2load,dir2use);
        D = load(['C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\Hindcast\Southern_LUT_Data\' fname]);
        [xg,yg] = findNearestGridPoint(D.X,D.Y,K(ii).lon,K(ii).lat); % Find closest wave grid point
        hsig(jj) = D.hs(xg,yg); % Grab wave point data
        spd = spd+2.5; %Increase Speed
    end
    x_axis = 5:2.5:5+(2.5*length(hsig)-1);
    ll(ii) = plot(x_axis,hsig,'o','MarkerSize',9,'MarkerFaceColor',b(ii,:),'MarkerEdgeColor',b(ii,:));
    hold on 
    line(x_axis,hsig,'Color',b(ii,:));
    hold on
    clear hsig 
    spd = 5;
end
grid on 
box on 
xlabel('Wind Speed [m/s]','FontSize',12)
ylabel('Hsig [m]','FontSize',12)
ax = gca;
ax.FontSize = 12;
lgd = legend(ll,arrayfun(@(x)(x.name),K,'un',0),...
    'interpreter','none',...
    'location','northwest');
lgd.FontSize = 12;

% -------------- Direction Sensitivity 
col = 1;
% direction = 'north';
b = parula(5);
dirWant = 0;
axes('position',[p_left+col*(p_wid+p_spacing)+.028 p_bot+row*(p_height+p_spacing) p_wid p_height])
for ii = 1:length(K)
    hsig = zeros(1,24);
    for jj = 1:length(hsig)
        spd2load = 15;
        dir2use = wrap2360(90-(dirWant - 180));
        fname = sprintf('southern_t200_s%d_d%d.mat',spd2load,dir2use);
        D = load(['C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\Hindcast\Southern_LUT_Data\' fname]);
        [xg,yg] = findNearestGridPoint(D.X,D.Y,K(ii).lon,K(ii).lat); % Find closest wave grid point
        hsig(jj) = D.hs(xg,yg); % Grab wave point data
        dirWant = dirWant + 15; %Increase Speed
    end
    x_axis = 0:15:345;
    ll(ii) = plot(x_axis,hsig,'o','MarkerSize',9,'MarkerFaceColor',b(ii,:),'MarkerEdgeColor',b(ii,:));
    hold on 
    line(x_axis,hsig,'Color',b(ii,:));
    hold on
    clear hsig 
    dirWant = 0;
end
grid on 
box on 
xlabel('Wind Direction [degrees]','FontSize',12)
ylabel(' ','FontSize',12)
ax = gca;
ax.FontSize = 12;
ax.XLim = [0 360];
% lgd = legend(ll,arrayfun(@(x)(x.name),K,'un',0),...
%     'interpreter','none',...
%     'location','north');
% lgd.FontSize = 12;


% ------------------ Water Level
col = 2;
axes('position',[p_left+col*(p_wid+p_spacing)+.06 p_bot+row*(p_height+p_spacing) p_wid p_height])
% direction = 'north';
switch direction   % 'Coming from the....' notion 
    case 'north'
        dirWant = 0;
        dir2use = wrap2360(90-(dirWant - 180));
    case 'south'
        dirWant = 180;
        dir2use = wrap2360(90-(dirWant - 180));
        fname1 = sprintf('southern_t00_s20_d%d.mat',dir2use);
        fname2 = 'SpatialFromSouth.mat';
    otherwise
        disp('Error with direction')
end
b = parula(5);
tide2load = -200;
for ii = 1:length(K)
    hsig = zeros(1,16);
    for jj = 1:length(hsig)
        spd2load = 15;    
        fname = sprintf('southern_t%i_s%d_d%d.mat',tide2load,spd2load,dir2use);
        if tide2load == 0
            tide2load = '00';
            fname = sprintf('southern_t%s_s%d_d%d.mat',tide2load,spd2load,dir2use);
            tide2load = 0;
        end
        D = load(['C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\Hindcast\Southern_LUT_Data\' fname]);
        [xg,yg] = findNearestGridPoint(D.X,D.Y,K(ii).lon,K(ii).lat); % Find closest wave grid point
        hsig(jj) = D.hs(xg,yg); % Grab wave point data
        tide2load = tide2load + 50; %Increase Speed
    end
    x_axis = -2:.5:5.5;
    ll(ii) = plot(x_axis,hsig,'o','MarkerSize',9,'MarkerFaceColor',b(ii,:),'MarkerEdgeColor',b(ii,:));
    hold on 
    line(x_axis,hsig,'Color',b(ii,:));
    hold on
    clear hsig 
    tide2load = -200;
end
grid on 
box on 
xlabel('Water Level [m]','FontSize',12)
ylabel(' ','FontSize',12)
ax = gca;
ax.FontSize = 12;
% lgd = legend(ll,arrayfun(@(x)(x.name),K,'un',0),...
%     'interpreter','none',...
%     'location','west');
% lgd.FontSize = 12;

fname = sprintf('WaveSensitivity_Speed_Direction_WaterLevel_WindsFromThe_%s_Speed_%d',upper(direction),spd2load);
printFig(gcf,fname,[12 11],'png',300)

%% Combined Figure of Direction and Wind Speed and Basemap of Locations 

% Specify Dimensions of Plot
p_left = .05;
p_right = .05;
p_top = .03;
p_bot = .07;
p_spacing = .01;
p_wid = ((1-p_right-p_left-p_spacing)/3);
p_height = (1-p_top-p_bot-p_spacing);

row = 0;
col = 0;
% -------------- Speed Sensitivity 
direction = 'north';
clf
switch direction   % 'Coming from the....' notion 
    case 'north'
        dirWant = 0;
        dir2use = wrap2360(90-(dirWant - 180));
    case 'south'
        dirWant = 180;
        dir2use = wrap2360(90-(dirWant - 180));
        fname1 = sprintf('southern_t00_s20_d%d.mat',dir2use);
        fname2 = 'SpatailFromSouth.mat';
    otherwise
        disp('Error with direction')
end
b = parula(5);
spd = 5;
axes('position',[p_left+col*(p_wid+p_spacing) p_bot+row*(p_height+p_spacing) p_wid p_height])
for ii = 1:length(K)
    hsig = zeros(1,10);
    for jj = 1:length(hsig)
        spd2load = round(spd);
        fname = sprintf('southern_t200_s%d_d%d.mat',spd2load,dir2use);
        D = load(['C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\Hindcast\Southern_LUT_Data\' fname]);
        [xg,yg] = findNearestGridPoint(D.X,D.Y,K(ii).lon,K(ii).lat); % Find closest wave grid point
        hsig(jj) = D.hs(xg,yg); % Grab wave point data
        spd = spd+2.5; %Increase Speed
    end
    x_axis = 5:2.5:5+(2.5*length(hsig)-1);
    ll(ii) = plot(x_axis,hsig,'o','MarkerSize',9,'MarkerFaceColor',b(ii,:),'MarkerEdgeColor',b(ii,:));
    hold on 
    line(x_axis,hsig,'Color',b(ii,:));
    hold on
    clear hsig 
    spd = 5;
end
grid on 
box on 
xlabel('Wind Speed [m/s]','FontSize',12)
ylabel('Hsig [m]','FontSize',12)
ax = gca;
ax.FontSize = 12;
lgd = legend(ll,arrayfun(@(x)(x.name),K,'un',0),...
    'interpreter','none',...
    'location','northwest');
lgd.FontSize = 12;

% -------------- Direction Sensitivity 
col = 1;
direction = 'south';
b = parula(5);
dirWant = 0;
axes('position',[p_left+col*(p_wid+p_spacing)+.028 p_bot+row*(p_height+p_spacing) p_wid p_height])
for ii = 1:length(K)
    hsig = zeros(1,24);
    for jj = 1:length(hsig)
        spd2load = 15;
        dir2use = wrap2360(90-(dirWant - 180));
        fname = sprintf('southern_t200_s%d_d%d.mat',spd2load,dir2use);
        D = load(['C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\Hindcast\Southern_LUT_Data\' fname]);
        [xg,yg] = findNearestGridPoint(D.X,D.Y,K(ii).lon,K(ii).lat); % Find closest wave grid point
        hsig(jj) = D.hs(xg,yg); % Grab wave point data
        dirWant = dirWant + 15; %Increase Speed
    end
    x_axis = 0:15:345;
    ll(ii) = plot(x_axis,hsig,'o','MarkerSize',9,'MarkerFaceColor',b(ii,:),'MarkerEdgeColor',b(ii,:));
    hold on 
    line(x_axis,hsig,'Color',b(ii,:));
    hold on
    clear hsig 
    dirWant = 0;
end
grid on 
box on 
xlabel('Wind Direction [degrees]','FontSize',12)
ylabel(' ','FontSize',12)
ax = gca;
ax.FontSize = 12;
% lgd = legend(ll,arrayfun(@(x)(x.name),K,'un',0),...
%     'interpreter','none',...
%     'location','north');
% lgd.FontSize = 12;


% ------------------ Basemap
col = 2;
axes('position',[p_left+col*(p_wid+p_spacing)+.028 p_bot+row*(p_height+p_spacing) p_wid p_height])
IM = load('C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\GeoTiffs\seattle_tacoma3.mat');
imagesc(IM.xm,IM.ym,IM.im)
set(gca,'ydir','normal')
hold on

for ii = 1:length(K)
    ll(ii)=plot(K(ii).x,K(ii).y,'o','MarkerSize',14,'MarkerFaceColor',b(ii,:),'MarkerEdgeColor','w','LineWidth',1.5);
end

xlabel('Easting [m]','FontSize',12)
ylabel('Northing [m]','FontSize',12)
ax = gca;
ax.FontSize = 12;

axis equal
xlim([5.2*10^5 5.55*10^5])
ylim([5.23*10^6 5.288*10^6])


lgd = legend(ll,arrayfun(@(x)(x.name),K,'un',0),...
    'interpreter','none',...
    'location','northwest');
lgd.FontSize = 12;








fname = sprintf('WaveSensitivty_Speed_Direction_SpeedOf%d_WithBasemap',spd2load);
printFig(gcf,fname,[11 8.5],'png',300)

% JUST BASEMAP
fname = sprintf('WaveSensitivty_Locations_Basemap',spd2load);
printFig(gcf,fname,[8.5 11],'png',300)

%% Other Way to Plot the Sensitivities - Stacked

% Specify Dimensions of Plot
p_left = .05;
p_right = .05;
p_top = .02;
p_bot = .055;
p_spacing = .055;
p_wid = ((1-p_right-p_left));
p_height = (1-p_top-p_bot-p_spacing)/3.2;

row = 2;
col = 0;
% -------------- Speed Sensitivity 
direction = 'north';
clf
switch direction   % 'Coming from the....' notion 
    case 'north'
        dirWant = 0;
        dir2use = wrap2360(90-(dirWant - 180));
    case 'south'
        dirWant = 180;
        dir2use = wrap2360(90-(dirWant - 180));
        fname1 = sprintf('southern_t00_s20_d%d.mat',dir2use);
        fname2 = 'SpatailFromSouth.mat';
    otherwise
        disp('Error with direction')
end
b = parula(5);
spd = 5;
axes('position',[p_left+col*(p_wid+p_spacing) p_bot+row*(p_height+p_spacing) p_wid p_height])
for ii = 1:length(K)
    hsig = zeros(1,10);
    for jj = 1:length(hsig)
        spd2load = round(spd);
        fname = sprintf('southern_t200_s%d_d%d.mat',spd2load,dir2use);
        D = load(['C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\Hindcast\Southern_LUT_Data\' fname]);
        [xg,yg] = findNearestGridPoint(D.X,D.Y,K(ii).lon,K(ii).lat); % Find closest wave grid point
        hsig(jj) = D.hs(xg,yg); % Grab wave point data
        spd = spd+2.5; %Increase Speed
    end
    x_axis = 5:2.5:5+(2.5*length(hsig)-1);
    ll(ii) = plot(x_axis,hsig,'o','MarkerSize',9,'MarkerFaceColor',b(ii,:),'MarkerEdgeColor',b(ii,:));
    hold on 
    line(x_axis,hsig,'Color',b(ii,:));
    hold on
    clear hsig 
    spd = 5;
end
grid on 
box on 
xlabel('Wind Speed [m/s]','FontSize',12)
ylabel('Hsig [m]','FontSize',12)
ax = gca;
ax.FontSize = 12;
ax.XLim = [5 28];
lgd = legend(ll,arrayfun(@(x)(x.name),K,'un',0),...
    'interpreter','none',...
    'location','northwest');
lgd.FontSize = 12;

% -------------- Direction Sensitivity 
row = 1;
% direction = 'north';
b = parula(5);
dirWant = 0;
axes('position',[p_left+col*(p_wid+p_spacing) p_bot+row*(p_height+p_spacing) p_wid p_height])
for ii = 1:length(K)
    hsig = zeros(1,24);
    for jj = 1:length(hsig)
        spd2load = 15;
        dir2use = wrap2360(90-(dirWant - 180));
        fname = sprintf('southern_t200_s%d_d%d.mat',spd2load,dir2use);
        D = load(['C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\Hindcast\Southern_LUT_Data\' fname]);
        [xg,yg] = findNearestGridPoint(D.X,D.Y,K(ii).lon,K(ii).lat); % Find closest wave grid point
        hsig(jj) = D.hs(xg,yg); % Grab wave point data
        dirWant = dirWant + 15; %Increase Speed
    end
    x_axis = 0:15:345;
    ll(ii) = plot(x_axis,hsig,'o','MarkerSize',9,'MarkerFaceColor',b(ii,:),'MarkerEdgeColor',b(ii,:));
    hold on 
    line(x_axis,hsig,'Color',b(ii,:));
    hold on
    clear hsig 
    dirWant = 0;
end
grid on 
box on 
xlabel('Wind Direction [degrees]','FontSize',12)
ylabel('Hsig [m]','FontSize',12)
ax = gca;
ax.FontSize = 12;
ax.XLim = [0 350];
% lgd = legend(ll,arrayfun(@(x)(x.name),K,'un',0),...
%     'interpreter','none',...
%     'location','north');
% lgd.FontSize = 12;


% ------------------ Water Level
row = 0;
axes('position',[p_left+col*(p_wid+p_spacing) p_bot+row*(p_height+p_spacing) p_wid p_height])
% direction = 'north';
switch direction   % 'Coming from the....' notion 
    case 'north'
        dirWant = 0;
        dir2use = wrap2360(90-(dirWant - 180));
    case 'south'
        dirWant = 180;
        dir2use = wrap2360(90-(dirWant - 180));
        fname1 = sprintf('southern_t00_s20_d%d.mat',dir2use);
        fname2 = 'SpatialFromSouth.mat';
    otherwise
        disp('Error with direction')
end
b = parula(5);
tide2load = -200;
for ii = 1:length(K)
    hsig = zeros(1,16);
    for jj = 1:length(hsig)
        spd2load = 15;    
        fname = sprintf('southern_t%i_s%d_d%d.mat',tide2load,spd2load,dir2use);
        if tide2load == 0
            tide2load = '00';
            fname = sprintf('southern_t%s_s%d_d%d.mat',tide2load,spd2load,dir2use);
            tide2load = 0;
        end
        D = load(['C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\Hindcast\Southern_LUT_Data\' fname]);
        [xg,yg] = findNearestGridPoint(D.X,D.Y,K(ii).lon,K(ii).lat); % Find closest wave grid point
        hsig(jj) = D.hs(xg,yg); % Grab wave point data
        tide2load = tide2load + 50; %Increase Speed
    end
    x_axis = -2:.5:5.5;
    ll(ii) = plot(x_axis,hsig,'o','MarkerSize',9,'MarkerFaceColor',b(ii,:),'MarkerEdgeColor',b(ii,:));
    hold on 
    line(x_axis,hsig,'Color',b(ii,:));
    hold on
    clear hsig 
    tide2load = -200;
end
grid on 
box on 
xlabel('Water Level [m]','FontSize',12)
ylabel('Hsig [m] ','FontSize',12)
ax = gca;
ax.FontSize = 12;
ax.XLim = [-2 5.6];
% lgd = legend(ll,arrayfun(@(x)(x.name),K,'un',0),...
%     'interpreter','none',...
%     'location','west');
% lgd.FontSize = 12;

fname = sprintf('WaveSensitivity_Speed_Direction_WaterLevel_WindsFromThe_%s_Speed_%d_STACKED',upper(direction),spd2load);
printFig(gcf,fname,[12 11],'png',300)