% Loads in MACA data and creates QQ Plots for specfic locations


% ------------------------ Load in MACA -----------------------------------
addpath C:\Functions_Matlab\time
addpath C:\Functions_Matlab
addpath C:\Functions_Matlab\cbrewer
addpath C:\Functions_Matlab\circ_stats
addpath C:\Functions_Matlab\mapping\kml
clearvars



temp_lat = ncread('E:\Abbas\Model_Met_Forcings\MACA\maca_u_GFDL_ESM2M_historical_1950_2005.nc','lat');
temp_lon = ncread('E:\Abbas\Model_Met_Forcings\MACA\maca_u_GFDL_ESM2M_historical_1950_2005.nc','lon');
temp_lon = -1*(360 - temp_lon);
[lat,lon] = meshgrid(temp_lat, temp_lon);


% ---------------- Load WA Coastline
load('C:\Users\ahooshmand\Desktop\PS_COSMOS\Salish_Model_Resources\WA_Spatial_Data\WA_coast.mat');
wa_lat = [];
wa_lon = [];
thresh = 200 ;
for j = 1:length(wa_coast)
    temp_x = wa_coast(j).X;
    temp_y = wa_coast(j).Y;
    if length(temp_x) >= thresh && j ~= 3348 % 3348 is oregon
        for m = 1:length(temp_x);
            wa_lat(end+1) = temp_y(m);
            wa_lon(end+1) = temp_x(m);
        end
    end
end

% --------------------------------------- Load in KML Points and find closest maca grid to it
kml_fol = 'C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\KML\PointLocations\';
k = dir('C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\KML\PointLocations\*.kml');
for ii = 1:length(k)
    K(ii) = kml2struct([kml_fol k(ii).name]);
    %     [K(ii).xg,K(ii).yg] = findNearestGridPoint(lon,lat,K(ii).Lon,K(ii).Lat); % Find closest wave grid point
end
% WHY THESE WON'T WORK IN THE SAME LOOP IS BEYOND ME!
for ii = 1:length(k)
    [K(ii).xg,K(ii).yg] = findNearestGridPoint(lon,lat,K(ii).Lon,K(ii).Lat); % Find closest wave grid point
end


% ---------------------- Historic Data ------------------------------------
Mh_u = ncread('E:\Abbas\Model_Met_Forcings\MACA\maca_u_GFDL_ESM2M_historical_1950_2005.nc','eastward_wind');
Mh_v = ncread('E:\Abbas\Model_Met_Forcings\MACA\maca_v_GFDL_ESM2M_historical_1950_2005.nc','northward_wind');
Mh_time = ncread('E:\Abbas\Model_Met_Forcings\MACA\maca_u_GFDL_ESM2M_historical_1950_2005.nc','time');
Mh_time = double(Mh_time); % convert to double for following step
Mh_time = datenum(1900,1,1+Mh_time); % Convert MACA time to Matlab Datenum
Mh_lat = ncread('E:\Abbas\Model_Met_Forcings\MACA\maca_u_GFDL_ESM2M_historical_1950_2005.nc','lat');
Mh_lon = ncread('E:\Abbas\Model_Met_Forcings\MACA\maca_u_GFDL_ESM2M_historical_1950_2005.nc','lon');
Mh_lon = -1*(360 - Mh_lon);
Mh_spd = hypot(Mh_u,Mh_v);
wnddir_temp = (180/pi)*atan2(Mh_u,Mh_v);
% Rotate winds to compass directions
wnddir_temp = 90 - wnddir_temp;
wnddir_temp(wnddir_temp<0)=wnddir_temp(wnddir_temp<0)+360;
% switch to conventional coming from dir
wnddir_temp = wnddir_temp+180;
Mh_dir = wrap2360(wnddir_temp);
clear wnddir_temp





% ---------------------- Future Data --------------------------------------
scenario = 'both';
switch scenario
    case 'rcp45'
        Mf_u = ncread('E:\Abbas\Model_Met_Forcings\MACA\maca_u_GFDL_ESM2M_rcp45_2006_2100.nc','eastward_wind');
        Mf_v = ncread('E:\Abbas\Model_Met_Forcings\MACA\maca_v_GFDL_ESM2M_rcp45_2006_2100.nc','northward_wind');
        Mf_time = ncread('E:\Abbas\Model_Met_Forcings\MACA\maca_u_GFDL_ESM2M_rcp45_2006_2100.nc','time');
        Mf_time = double(Mf_time); % convert to double for following step
        Mf_time = datenum(1900,1,1+Mf_time); % Convert MACA time to Matlab Datenum
        Mf_spd = hypot(Mf_u,Mf_v);
        wnddir_temp = (180/pi)*atan2(Mf_u,Mf_v);
        % Rotate winds to compass directions
        wnddir_temp = 90 - wnddir_temp;
        wnddir_temp(wnddir_temp<0)=wnddir_temp(wnddir_temp<0)+360;
        % switch to conventional coming from dir
        wnddir_temp = wnddir_temp+180;
        Mf_dir = wrap2360(wnddir_temp);
        clear wnddir_temp
        
        % Concatenate Matrices
        M.u = cat(3,Mh_u,Mf_u);
        clear Mh_u Mf_u
        M.v = cat(3,Mh_v,Mf_v);
        clear Mh_v Mf_v
        M.time = cat(1,Mh_time,Mf_time);
        clear Mh_time Mf_time
        [M.lat,M.lon] = meshgrid(Mh_lat, Mh_lon);
        clear Mh_lat Mh_lon
        M.spd = cat(3,Mh_spd,Mf_spd);
        clear Mh_spd Mf_spd
        M.dir = cat(3,Mh_dir,Mf_dir);
        clear Mh_dir Mf_dir
    case 'rcp85'
        Mf_u = ncread('E:\Abbas\Model_Met_Forcings\MACA\maca_u_GFDL_ESM2M_rcp85_2006_2100.nc','eastward_wind');
        Mf_v = ncread('E:\Abbas\Model_Met_Forcings\MACA\maca_v_GFDL_ESM2M_rcp85_2006_2100.nc','northward_wind');
        Mf_time = ncread('E:\Abbas\Model_Met_Forcings\MACA\maca_u_GFDL_ESM2M_rcp85_2006_2100.nc','time');
        Mf_time = double(Mf_time); % convert to double for following step
        Mf_time = datenum(1900,1,1+Mf_time); % Convert MACA time to Matlab Datenum
        Mf_spd = hypot(Mf_u,Mf_v);
        wnddir_temp = (180/pi)*atan2(Mf_u,Mf_v);
        % Rotate winds to compass directions
        wnddir_temp = 90 - wnddir_temp;
        wnddir_temp(wnddir_temp<0)=wnddir_temp(wnddir_temp<0)+360;
        % switch to conventional coming from dir
        wnddir_temp = wnddir_temp+180;
        Mf_dir = wrap2360(wnddir_temp);
        clear wnddir_temp
        
        % Concatenate Matrices
        M.u = cat(3,Mh_u,Mf_u);
        clear Mh_u Mf_u
        M.v = cat(3,Mh_v,Mf_v);
        clear Mh_v Mf_v
        M.time = cat(1,Mh_time,Mf_time);
        clear Mh_time Mf_time
        [M.lat,M.lon] = meshgrid(Mh_lat, Mh_lon);
        clear Mh_lat Mh_lon
        M.spd = cat(3,Mh_spd,Mf_spd);
        clear Mh_spd Mf_spd
        M.dir = cat(3,Mh_dir,Mf_dir);
        clear Mh_dir Mf_dir
    case 'both'
        Mf_u = ncread('E:\Abbas\Model_Met_Forcings\MACA\maca_u_GFDL_ESM2M_rcp45_2006_2100.nc','eastward_wind');
        Mf_v = ncread('E:\Abbas\Model_Met_Forcings\MACA\maca_v_GFDL_ESM2M_rcp45_2006_2100.nc','northward_wind');
        Mf_time = ncread('E:\Abbas\Model_Met_Forcings\MACA\maca_u_GFDL_ESM2M_rcp45_2006_2100.nc','time');
        Mf_time = double(Mf_time); % convert to double for following step
        Mf_time = datenum(1900,1,1+Mf_time); % Convert MACA time to Matlab Datenum
        Mf_spd = hypot(Mf_u,Mf_v);
        wnddir_temp = (180/pi)*atan2(Mf_u,Mf_v);
        % Rotate winds to compass directions
        wnddir_temp = 90 - wnddir_temp;
        wnddir_temp(wnddir_temp<0)=wnddir_temp(wnddir_temp<0)+360;
        % switch to conventional coming from dir
        wnddir_temp = wnddir_temp+180;
        Mf_dir = wrap2360(wnddir_temp);
        clear wnddir_temp
        
        % Concatenate Matrices
        M4.u = cat(3,Mh_u,Mf_u);
        clear Mf_u
        M4.v = cat(3,Mh_v,Mf_v);
        clear Mf_v
        M4.time = cat(1,Mh_time,Mf_time);
        clear Mf_time
        [M4.lat,M4.lon] = meshgrid(Mh_lat, Mh_lon);
        M4.spd = cat(3,Mh_spd,Mf_spd);
        clear Mf_spd
        M4.dir = cat(3,Mh_dir,Mf_dir);
        clear Mf_dir
        
        for ii = 1:length(K)
            K(ii).u4 = squeeze(M4.u(K(ii).xg,K(ii).yg,:));
            K(ii).v4 = squeeze(M4.v(K(ii).xg,K(ii).yg,:));
            K(ii).spd4 = squeeze(M4.spd(K(ii).xg,K(ii).yg,:));
            K(ii).dir4 = squeeze(M4.dir(K(ii).xg,K(ii).yg,:));
            K(ii).time = M4.time;
            
        end
        clear M4
        
        % Now 8.5
        Mf_u = ncread('E:\Abbas\Model_Met_Forcings\MACA\maca_u_GFDL_ESM2M_rcp85_2006_2100.nc','eastward_wind');
        Mf_v = ncread('E:\Abbas\Model_Met_Forcings\MACA\maca_v_GFDL_ESM2M_rcp85_2006_2100.nc','northward_wind');
        Mf_time = ncread('E:\Abbas\Model_Met_Forcings\MACA\maca_u_GFDL_ESM2M_rcp85_2006_2100.nc','time');
        Mf_time = double(Mf_time); % convert to double for following step
        Mf_time = datenum(1900,1,1+Mf_time); % Convert MACA time to Matlab Datenum
        Mf_spd = hypot(Mf_u,Mf_v);
        wnddir_temp = (180/pi)*atan2(Mf_u,Mf_v);
        % Rotate winds to compass directions
        wnddir_temp = 90 - wnddir_temp;
        wnddir_temp(wnddir_temp<0)=wnddir_temp(wnddir_temp<0)+360;
        % switch to conventional coming from dir
        wnddir_temp = wnddir_temp+180;
        Mf_dir = wrap2360(wnddir_temp);
        clear wnddir_temp
        
        % Concatenate Matrices
        M8.u = cat(3,Mh_u,Mf_u);
        clear Mh_u Mf_u
        M8.v = cat(3,Mh_v,Mf_v);
        clear Mh_v Mf_v
        M8.time = cat(1,Mh_time,Mf_time);
        clear Mh_time Mf_time
        [M8.lat,M8.lon] = meshgrid(Mh_lat, Mh_lon);
        clear Mh_lat Mh_lon
        M8.spd = cat(3,Mh_spd,Mf_spd);
        clear Mh_spd Mf_spd
        M8.dir = cat(3,Mh_dir,Mf_dir);
        clear Mh_dir Mf_dir
        
        for ii = 1:length(K)
            K(ii).u8 = squeeze(M8.u(K(ii).xg,K(ii).yg,:));
            K(ii).v8 = squeeze(M8.v(K(ii).xg,K(ii).yg,:));
            K(ii).spd8 = squeeze(M8.spd(K(ii).xg,K(ii).yg,:));
            K(ii).dir8 = squeeze(M8.dir(K(ii).xg,K(ii).yg,:));
            
        end
        clear M8
        
end

Mh_u = ncread('E:\Abbas\Model_Met_Forcings\MACA\maca_u_GFDL_ESM2M_historical_1950_2005.nc','eastward_wind'); % THIS IS FOR PLOTTING LATER
return
%% QQ Plot for MACA

hist = find(year(K(1).time) >= 1950 & year(K(1).time) <= 2020);
fut = find(year(K(1).time) >= 2030 & year(K(1).time) <= 2100);
hist(1:50) = [];

scenario = 'rcp85'
clf
for ii = 1:length(K)
    
    switch scenario
        case 'rcp45'
            data_hist = K(ii).spd4(hist);
            data_fut = K(ii).spd4(fut);
            fname = sprintf('QQ_%s_RCP45',K(ii).Name);
        case 'rcp85'
            data_hist = K(ii).spd8(hist);
            data_fut = K(ii).spd8(fut);
            fname = sprintf('QQ_%s_RCP85',K(ii).Name);
    end
    
    cdf1 = sort(data_hist,'ascend');
    cdf1y = linspace(0,1,length(cdf1));
    
    cdf2 = sort(data_fut,'ascend');
    cdf2y = linspace(0,1,length(cdf2));
    
    %Interp cdfs
%     cdf2i = interp1(cdf2y,cdf2,cdf1y);
    
    
    x_axis = linspace(0,1,length(cdf1));
    clf
    plot(cdf1,cdf2,'r+')
    tit = sprintf('%s',K(ii).Name);
    title(tit)
    grid on
    hold on
    line([0 50],[0 50],'Color','k','LineStyle','--')
    xlabel('MACA Historic: 1950 - 2020');
    ylabel('MACA Future: 2030 - 2100');
    xlim([0 18])
    ylim([0 18])
    
    set(gca,'FontSize',14)
    printFig(gcf,fname,[14 14],'png',300)
end

% printFig(gcf,'MACA_85_QQ_Speed',[14 14],'png',300)

%% Plot WA Basemap and Location of each Grid Cell and KML Point
clf
plot(wa_lon,wa_lat)
hold on
b = parula(length(K));
for ii = 1:length(K)
    c = plot(K(ii).Lon,K(ii).Lat,'o','MarkerSize',8,'MarkerFaceColor','w','MarkerEdgeColor','k');
    s = plot(temp_lon(K(ii).xg),temp_lat(K(ii).yg),'r*','MarkerSize',8);
end
xlim([-125 -122])
ylim([46.5 49])
xlabel('Degrees Longitude','FontSize',14)
ylabel('Degrees Latitude','FontSize',14)
set(gca,'fontsize',14)

grid on
g1 = plot(lon,lat,'.','Color',[.7 .7 .7],'MarkerSize',3);
temp = Mh_u(:,:,1);
temp = isnan(temp);
g2 = plot(lon(~isnan(Mh_u(:,:,1))),lat(~isnan(Mh_u(:,:,1))),'k.','MarkerSize',4);

legend([c s g2],'KML Point','Closest MACA Point', 'MACA Data Grid','location','northwest');
printFig(gcf,'MACA_Grid_KMLPoints',[14 14],'png',300)
%% QQ Plot of Differing Time Values and Histogram

hist = find(year(K(1).time) >= 1950 & year(K(1).time) <= 2020);
fut = find(year(K(1).time) >= 2020 & year(K(1).time) <= 2100);

scenario = 'rcp45';

for ii = 1:length(K)
    switch scenario 
        case 'rcp45'
            data_hist = K(ii).spd4(hist);
            data_fut = K(ii).spd4(fut);
            fname = sprintf('QQ_Hist_%s_RCP45',K(ii).Name);
        case 'rcp85'
            data_hist = K(ii).spd8(hist);
            data_fut = K(ii).spd8(fut);
            fname = sprintf('QQ_Hist_%s_RCP85',K(ii).Name);
    end
    
    cdf1 = sort(data_hist,'ascend');
    cdf1y = linspace(0,1,length(cdf1));
    
    cdf2 = sort(data_fut,'ascend');
    cdf2y = linspace(0,1,length(cdf2));
    
    %Interp cdfs
    cdf2i = interp1(cdf2y,cdf2,cdf1y);
    
    clf
    edges = 0:.25:15;
    bins = conv(edges,.5*[1 1],'valid');
    n1 = histcounts(cdf1,edges,'Normalization','probability');
    n2 = histcounts(cdf2,edges,'Normalization','probability');
    subplot(121)
    plot(bins,n1,bins,n2)
    ylabel('Probability','FontSize',14)
    xlabel('Wind Speed [m/s]','FontSize',14)
    set(gca,'fontsize',14)
    legend('1950-2020','2020-2100')

    
    subplot(122)
    plot(cdf1,cdf2i)
    grid on
    hold on
    plot([0 50],[0 50],'-k')
    box on
    xlabel('Speed [m/s]: 1950-2020','FontSize',14)
    ylabel('Speed [m/s]: 2020-2100','FontSize',14)
    set(gca,'fontsize',14)
    xlim([0 18])
    ylim([0 18])
    
    printFig(gcf,fname,[14 14],'png',300)
end

%% QQ Plots of 4.5 vs 8.5


fut = find(year(K(1).time) >= 2020 & year(K(1).time) <= 2100);

scenario = 'rcp45';

for ii = 1:length(K)


    data_45 = K(ii).spd4(fut);
    
    data_85 = K(ii).spd8(fut);

    cdf1 = sort(data_45,'ascend'); % 4.5
    cdf1y = linspace(0,1,length(cdf1));
    
    cdf2 = sort(data_85,'ascend'); % 8.5
    cdf2y = linspace(0,1,length(cdf2));
    
    %Interp cdfs
%     cdf2i = interp1(cdf2y,cdf2,cdf1y);
    
    clf
    edges = 0:.25:15;
    bins = conv(edges,.5*[1 1],'valid');
    n1 = histcounts(cdf1,edges,'Normalization','probability');
    n2 = histcounts(cdf2,edges,'Normalization','probability');
    subplot(121)
    plot(bins,n1,bins,n2)
    ylabel('Probability','FontSize',14)
    xlabel('Wind Speed [m/s]','FontSize',14)
    set(gca,'fontsize',14)
    legend('1950-2020','2020-2100')

    
    subplot(122)
    plot(cdf1,cdf2,'r+') 
    grid on
    hold on
    plot([0 50],[0 50],'-k')
    box on
    xlabel('Speed [m/s]: RCP 4.5','FontSize',14)
    ylabel('Speed [m/s]: RCP 8.5','FontSize',14)
    set(gca,'fontsize',14)
    xlim([0 18])
    ylim([0 18])
    fname = sprintf('%s_45_vs_85',K(ii).Name);
    printFig(gcf,fname,[14 14],'png',300)
end

%% MultiPannel QQ Plots 

p_left = .05;
p_right = .05;
p_top = .03;
p_bot = .07;
p_spacing = .01;
p_wid = ((1-p_right-p_left-p_spacing)/3);
p_height = (1-p_top-p_bot-p_spacing)/3;

row = 2;
col = 0;

% hist = find(year(K(1).time) >= 1950 & year(K(1).time) <= 2020);
% fut = find(year(K(1).time) >= 2020 & year(K(1).time) <= 2100);
% hist(1:50) = [];

scenario = '85Maxs';
clf
for ii = 1:length(K)
    
    switch scenario
        case 'rcp45'
            data_hist = K(ii).spd4(hist);
            data_fut = K(ii).spd4(fut);
            fname = sprintf('QQ_%s_RCP45',K(ii).Name);
        case 'rcp85'
            data_hist = K(ii).spd8(hist);
            data_fut = K(ii).spd8(fut);
            fname = sprintf('QQ_%s_RCP85',K(ii).Name);
        case '45vs85'
            data_hist = K(ii).spd4(fut);
            data_fut = K(ii).spd8(fut);
        case '85Maxs'
            hist = 1:71;
            fut = 51:151;
            data_hist = squeeze(spd_max(K(ii).xg,K(ii).yg,hist));
            data_fut = squeeze(spd_max(K(ii).xg,K(ii).yg,fut));
            fname = sprintf('QQ_Hist_%s_RCP85_MAXSPEEDS',K(ii).Name);
    end
    
    cdf1 = sort(data_hist,'ascend');
    cdf1y = linspace(0,1,length(cdf1));
    
    cdf2 = sort(data_fut,'ascend');
    cdf2y = linspace(0,1,length(cdf2));
    
    %Interp cdfs
    cdf2i = interp1(cdf2y,cdf2,cdf1y);
    
    
    x_axis = linspace(0,1,length(cdf1));
    axes('position',[p_left+col*(p_wid+p_spacing) p_bot+row*(p_height+p_spacing) p_wid p_height])
    plot(cdf1,cdf2i,'r+')
    tit = sprintf('%s',K(ii).Name);
    title(tit)
    set(get(gca,'title'),'Position',[7 0.4 1.00011])
    grid on
    hold on
    line([0 50],[0 50],'Color','k','LineStyle','--')
    if ii == 2 || ii == 3 || ii == 5 || ii == 6
        xlabel(' ')
        ylabel(' ')
        set(gca,'XTickLabel',[])
        set(gca,'YTickLabel',[])
        xlim([0 18])
        ylim([0 18])
        switch scenario
            case '85Maxs'
                xlim([5 17])
                ylim([5 17])
                set(get(gca,'title'),'Position',[11 6.5 1.00011])
        end
    elseif ii == 1 || ii == 4
        switch scenario
            case 'rcp45'
                xlabel(' ')
                ylabel('MACA 2020 - 2100 [m/s]')
            case 'rcp85'
                xlabel(' ')
                ylabel('MACA 2020 - 2100 [m/s]')
            case '45vs85'
                xlabel(' ')
                ylabel('RCP 8.5 2020-2100 [m/s]')
            case '85Maxs'
                xlabel(' ')
                ylabel('Max Spd. 2020-2100 [m/s]')
        end
        set(gca,'XtickLabel',[])
        xlim([0 18])
        ylim([0 18])
        switch scenario
            case '85Maxs'
                xlim([5 17])
                ylim([5 17])
                set(get(gca,'title'),'Position',[11 6.5 1.00011])
        end
    elseif ii == 8 || ii == 9
        switch scenario 
            case 'rcp45'
                ylabel(' ')
                xlabel('MACA 1950 - 2020 [m/s]')
            case 'rcp85'
                ylabel(' ')
                xlabel('MACA 1950 - 2020 [m/s]')
            case '45vs85'
                ylabel(' ')
                xlabel('RCP 4.5 2020-2100 [m/s]')
            case '85Maxs'
                ylabel(' ')
                xlabel('Max Spd. 1950-2020 [m/s]')
        end
        set(gca,'YtickLabel',[])
        xlim([0 18])
        ylim([0 18])
        switch scenario
            case '85Maxs'
                xlim([5 17])
                ylim([5 17])
                set(get(gca,'title'),'Position',[11 6.5 1.00011])
        end
    else
        switch scenario
            case 'rcp45'
                xlabel('MACA 1950 - 2020 [m/s]');
                ylabel('MACA 2020 - 2100 [m/s]');
            case 'rcp85'
                xlabel('MACA 1950 - 2020 [m/s]');
                ylabel('MACA 2020 - 2100 [m/s]');
            case '45vs85'
                xlabel('RCP 4.5 2020-2100 [m/s]')
                ylabel('RCP 8.5 2020-2100 [m/s]')
            case '85Maxs'
               xlabel('Max Spd. 1950-2020 [m/s]')
               ylabel('Max Spd. 2020-2100 [m/s]')
        end   
        xlim([0 18])
        ylim([0 18])
        switch scenario
            case '85Maxs'
                xlim([5 17])
                ylim([5 17])
                set(get(gca,'title'),'Position',[11 6.5 1.00011])
        end
    end
    set(gca,'FontSize',14)
    
    
    % Increment position 
    col = col + 1;
    if col >=3
        col = 0;
        row = row - 1;
    end
end
fname = sprintf('QQ_Multipanel_%s',scenario);
printFig(gcf,fname,[14 14],'png',300)

%% QQ Plot of Max Speeds for varying times

load('MACA_Stats_Data_rcp85_MaxsGrabbed_1.mat');
years = 1950:1:2100;
ind = find(years == 2020);
hist = 1:ind;
fut = ind:length(years);
for ii = 2:length(K)

    data_hist = squeeze(spd_max(K(ii).xg,K(ii).yg,hist));
    data_fut = squeeze(spd_max(K(ii).xg,K(ii).yg,fut));
    fname = sprintf('QQ_Hist_%s_RCP85_MAXSPEEDS',K(ii).Name);

    [cdf1,Ih] = sort(data_hist,'ascend');
    cdf1y = linspace(0,1,length(cdf1));
    
    [cdf2,If] = sort(data_fut,'ascend');
    cdf2y = linspace(0,1,length(cdf2));
    
    %Interp cdfs
    cdf2i = interp1(cdf2y,cdf2,cdf1y);
    
    clf
    edges = 0:.25:15;
    bins = conv(edges,.5*[1 1],'valid');
    n1 = histcounts(cdf1,edges,'Normalization','probability');
    n2 = histcounts(cdf2,edges,'Normalization','probability');
    subplot(121)
    plot(bins,n1,bins,n2)
    ylabel('Probability','FontSize',14)
    xlabel('Wind Speed [m/s]','FontSize',14)
    set(gca,'fontsize',14)
    legend('1950-2020','2020-2100')
    xlim([7 16])

    
    subplot(122)
    plot(cdf1,cdf2i,'r+')
    grid on
    hold on
    plot([0 50],[0 50],'-k')
    box on
    xlabel('Yearly Mean Max Speed [m/s]: 1950-2020','FontSize',14)
    ylabel('Yearly Mean Max Speed [m/s]: 2020-2100','FontSize',14)
    set(gca,'fontsize',14)
    xlim([5 18])
    ylim([5 18])
    
%     printFig(gcf,fname,[14 14],'png',300)
pause
end
