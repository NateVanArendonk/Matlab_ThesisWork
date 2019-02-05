% What this function does 
% 1. Plots Climate windoes for MACA
% 2. Plots 2 long period climate windows for MACA
% 3. Plots of percent differences between time periods 
% 4. Compare Distributions
% 5. Plots The average wind speed for climate period and performs
%       regression
% 6. Plots of Binned wind speeds and directions 
% 7. Heat Map of Speed through time
% 8. Directional differences over time
% 9. Mean of 1950 - present compared to future projections 
% 10. 90th percentile of winds and how they change in the future
% 11. QQ Plots for MACA 


% IMPORTANT READ WHAT SCRIPT DOES ABOVE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



% ------------------------ Load in MACA -----------------------------------
addpath C:\Functions_Matlab\time
addpath C:\Functions_Matlab
addpath C:\Functions_Matlab\cbrewer
addpath C:\Functions_Matlab\circ_stats
clearvars
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
scenario = 'rcp85';
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
    case 'rcp85'
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
end

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

% Load in WA Coast 
load('WA_coast.mat');
% Find all shapes above threshold
wa_lat = [];
wa_lon = [];
thresh = 500 ;
for j = 1:length(wa_coast)
    temp_x = wa_coast(j).X;
    temp_y = wa_coast(j).Y;
    if length(temp_x) >= thresh %&& j ~= 3348 % 3348 is oregon
        for m = 1:length(temp_x);
            wa_lat(end+1) = temp_y(m);
            wa_lon(end+1) = temp_x(m);
        end
    end
end
clear j m temp_x temp_y thresh wa_coast
return
%% Grab MACA time segments - yearly vs long period 

% Time Periods
decade_starts = 1950:10:2070;
decade_ends = 1980:10:2100;
years = 1950:1:2100;

for xx = 1:length(decade_starts)
    % Find MACA time in each decade and years from year variable in each decade
    decade{xx} = find(year(M.time) >= decade_starts(xx) & year(M.time) <= decade_ends(xx));
    dec_yrs{xx} = find(years >= decade_starts(xx) & years <= decade_ends(xx));
end

% Find Time indices of variables 
period1 = find(year(M.time) >= 1950 & year(M.time) <= 2020);
period2 = find(year(M.time) >= 2020 & year(M.time) <= 2100);
years1 = find(years >= 1950 & years <= 2020);
years2 = find(years >= 2020 & years <= 2100);
years3 = find(years >= 2030 & years <= 2100);

% Allocate
[r,c] = size(M.u(:,:,1));
spd_max = zeros(r,c,length(years)); 
spd_avg = spd_max;
u_avg = spd_max;
v_avg = spd_max;
u_max = spd_max;
v_max = spd_max;
dir_max_mean = spd_max; % Pre-Allocate
dir_avg = spd_max;

max2grab = 5; 
% number of maxes to grab for average max per year
% Grab n-max values per year for each grid cell, and then average them to
% represent the 'max' per year per grid cell.  Also calculate the average
% direction for the maxes.  

for y = 1:length(years)
    temp_spd = M.spd(:,:,year(M.time)==years(y)); % Get all the speed indices that are for that year
    temp_dir = M.dir(:,:,year(M.time)==years(y)); % Get all the direction indices that are for that year
    temp_u = M.u(:,:,year(M.time)==years(y)); % Get all of U directions for that year
    temp_v = M.v(:,:,year(M.time)==years(y)); % Get all of V directions for that year
    for m = 1:r % Go row by row
        for n = 1:c % Go column by column 
            spds = temp_spd(m,n,:); spds = spds(:);% Grab the specific speed at that row and column for all times in the year and vectorize
            dirs = temp_dir(m,n,:); dirs = dirs(:);% Grab the specific directions at that row and column for all times in the year and vecotirze 
            uu = temp_u(m,n,:); uu = uu(:);% Grab the specific U directions at that row and column for all times in the year and vecotirze 
            vv = temp_v(m,n,:); vv = vv(:);% Grab the specific V directions at that row and column for all times in the year and vecotirze 
            nan_inds = find(isnan(spds)); % Find all the indices that are NaNs
            if length(nan_inds) ~= 365 % if the location isn't a NaN's only spot
                [B,I] = sort(spds,'descend'); % Sort the speeds with highest at the top
                spd_max(m,n,y) = mean(B(1:max2grab)); % Populate the matrix with the average of the top 'n' speeds
                dir_max_mean(m,n,y) = wrap2360(rad2deg(circ_mean(deg2rad(dirs(I(1:max2grab))),[],1))); % Take the average of the top N directions
                dir_avg(m,n,y) = wrap2360(rad2deg(circ_mean(deg2rad(dirs),[],1))); % Take the average of all directions for that year
                u_max(m,n,y) = mean(uu(I(1:max2grab))); % Take the average of the N max U values
                v_max(m,n,y) = mean(vv(I(1:max2grab))); % Take the average of the N max V values 
                spd_avg(m,n,y) = mean(spds);
                u_avg(m,n,y) = mean(uu);
                v_avg(m,n,y) = mean(vv);
                spd_90(m,n,y) = prctile(spds,90);
                spd_95(m,n,y) = prctile(spds,95);
                spd_99(m,n,y) = prctile(spds,99);
            else % If the point is all NaNs
                spd_max(m,n,y) = NaN;
                dir_max_mean(m,n,y) = NaN;
                u_max(m,n,y) = NaN;
                v_max(m,n,y) = NaN;
                dir_avg(m,n,y) = NaN;
                spd_avg(m,n,y) = NaN;
                u_avg(m,n,y) = NaN;
                v_avg(m,n,y) = NaN;
                spd_90(m,n,y) = NaN;
                spd_95(m,n,y) = NaN;
                spd_99(m,n,y) = NaN;
            end
        end
    end
end
clear c m n r uu vv xx yy 
fname = sprintf('MACA_Stats_Data_%s',scenario');
lat = M.lat; lon = M.lon;
save(fname,'spd_max','dir_max_mean','dir_avg','u_max','v_max','u_avg','v_avg','spd_90','spd_95','spd_99','lat','lon');
load('MACA_Stats_Data_rcp85.mat');
return
%% Compare MACA over Long Time Periods by Distribution
clf
% Plot distributions of averages - Non_max speed averages
xbins = 0:.5:7;
clf
histogram(spd_avg(:,:,years1),xbins,'Normalization','probability')
hold on
histogram(spd_avg(:,:,years2),xbins,'Normalization','probability','FaceAlpha',0.4)
set(gca,'FontSize',18)
xlabel('Average Wind Speed [m/s]')
ylabel('Probability')
legend('1950 - 2020','2020 - 2100')
printFig(gcf, 'MACA_Compare_Distributions_MEAN', [14 14], 'png', 150)

% Plot distributions of Maxes
xbins = 5:1:25;
clf
histogram(nan_mean(spd_max(:,:,years1),3),xbins,'Normalization','probability')
hold on
histogram(nan_mean(spd_max(:,:,years2),3),xbins,'Normalization','probability','FaceAlpha',0.4)
set(gca,'FontSize',18)
xlabel('Max Wind Speed [m/s]')
ylabel('Probability')
legend('1950 - 2020','2020 - 2100')
printFig(gcf, 'MACA_Compare_Distributions_MAX', [14 14], 'png', 150)

%% Compare MACA over Long Time Periods (1950 - 2020 & 2020 - 2100) with Pcolor

left_bound = min(M.lon(:));
right_bound = -122.1;
metric = 'max';
switch metric
    case 'mean'
        % First Plot
        clf
        pcolor(M.lon,M.lat,nan_mean(spd_avg(:,:,years1),3));
        axis equal
        xlim([left_bound right_bound])
        shading interp
        hold on
        % Add WA Coastline
        plot(wa_lon, wa_lat, 'Color', [1 1 1])
        hold on
        % Add Quiver
        temp_u = nan_mean(u_avg(:,:,years1),3);
        temp_v = nan_mean(v_avg(:,:,years1),3);
        [r,c] = size(temp_u);
        sub = 1:29:r*c;
        quiver(M.lon(sub), M.lat(sub), temp_u(sub),temp_v(sub),'Color','k');
        title('MACA: 1950 - 2020 Mean Windspeed and Direction')
        xlabel('Degrees Longitude')
        ylabel('Degrees Latitude')
        chan = colorbar('Location','EastOutSide');
        caxis([0 6])
        set(gca,'FontSize',18)
        printFig(gcf, 'MACA_1950_2020_Compare_MEAN', [14 14], 'png', 150)
        %pause
        
        % Second Plot
        clf
        pcolor(M.lon,M.lat,nan_mean(spd_avg(:,:,years2),3));
        axis equal
        xlim([left_bound right_bound])
        shading interp
        hold on
        % Add WA Coastline
        plot(wa_lon, wa_lat, 'Color', [1 1 1])
        hold on
        % Add Quiver
        temp_u = nan_mean(u_avg(:,:,years2),3);
        temp_v = nan_mean(v_avg(:,:,years2),3);
        [r,c] = size(temp_u);
        sub = 1:29:r*c;
        quiver(M.lon(sub), M.lat(sub), temp_u(sub),temp_v(sub),'Color','k');
        xlabel('Degrees Longitude')
        ylabel('Degrees Latitude')
        chan = colorbar('Location','EastOutSide');
        caxis([0 6])
        set(gca,'FontSize',18)
        title('MACA: 2020 - 2100 Mean Windspeed and Direction')
        printFig(gcf, 'MACA_2020_2100_Compare_MEAN', [14 14], 'png', 150)
    
    case 'max' 
        % First Plot
        clf
        pcolor(M.lon,M.lat,nan_mean(spd_max(:,:,years1),3));
        axis equal
        xlim([left_bound right_bound])
        shading interp
        hold on
        % Add WA Coastline
        plot(wa_lon, wa_lat, 'Color', [1 1 1])
        hold on
        % Add Quiver
        temp_u = nan_mean(v_max(:,:,years1),3);
        temp_v = nan_mean(v_max(:,:,years1),3);
        [r,c] = size(temp_u);
        sub = 1:29:r*c;
        quiver(M.lon(sub), M.lat(sub), temp_u(sub),temp_v(sub),'Color','k');
        title('MACA: 1950 - 2020 Max Wind speed and Direction')
        xlabel('Degrees Longitude')
        ylabel('Degrees Latitude')
        chan = colorbar('Location','EastOutSide');
        caxis([0 18])
        set(gca,'FontSize',18)
        title('MACA: 1950 - 2020 Max Windspeed and Direction')
        printFig(gcf, 'MACA_1950_2020_Compare_MAX', [14 14], 'png', 150)
        
        % Second time period Plot
        clf
        pcolor(M.lon,M.lat,nan_mean(spd_max(:,:,years2),3));
        axis equal
        xlim([left_bound right_bound])
        shading interp
        hold on
        % Add WA Coastline
        plot(wa_lon, wa_lat, 'Color', [1 1 1])
        hold on
        % Add Quiver
        temp_u = nan_mean(v_max(:,:,years2),3);
        temp_v = nan_mean(v_max(:,:,years2),3);
        [r,c] = size(temp_u);
        sub = 1:29:r*c;
        quiver(M.lon(sub), M.lat(sub), temp_u(sub),temp_v(sub),'Color','k');
        hold on
        title('MACA: 2020 - 2100 Max Wind speed and Direction')
        xlabel('Degrees Longitude')
        ylabel('Degrees Latitude')
        chan = colorbar('Location','EastOutSide');
        caxis([0 20])
        set(gca,'FontSize',18)
        title('MACA: 2020 - 2100 Max Windspeed and Direction')
        printFig(gcf, 'MACA_2020_2100_Compare_MAX', [14 14], 'png', 150)
end
%% Plot Percent difference between the time periods 
metric = 'max';
param = 'absolute';
% First between the two longer periods then between the decadal periods 

% disp('Max')
% disp([min(var_plot_max(:)) max(var_plot_max(:))])
% disp('Mean')
% disp([min(var_plot_mean(:)) max(var_plot_mean(:))])

clf
N_c = 100;
mycolors = flipud(cbrewer('div','RdBu',N_c));
colormap(mycolors)

switch metric 
    case 'max'
        switch param
            case 'relative'
                var_plot_max = 100*(((nan_mean(spd_max(:,:,years2),3)-nan_mean(spd_max(:,:,years1),3))./nan_mean(spd_max(:,:,years2),3)));
                pcolor(M.lon,M.lat,var_plot_max)
                shading flat
                c = colorbar;
                c.Label.String = 'Percent Difference - Max Wind Speed [m/s]';
                caxis([-5 5]);
            case 'absolute'
                var_plot_max = nan_mean(spd_max(:,:,years2),3)-nan_mean(spd_max(:,:,years1),3);
%                 var_plot_max = 100*(((nan_mean(spd_max(:,:,years1),3)-nan_mean(spd_max(:,:,years2),3))./nan_mean(spd_max(:,:,years1),3)));
                pcolor(M.lon,M.lat,var_plot_max)
                shading flat
                c = colorbar;
                c.Label.String = 'Difference in Max Wind Speed [m/s]';
                caxis([-1 1]);
        end
        title('Difference in Max Wind Speeds: 1950 - 2020 vs 2020 - 2100')
        hold on
        plot(wa_lon, wa_lat,'Color','k')
        axis equal
        xlim([min(M.lon(:)) -122.1])
        xlabel('Degrees Longitude')
        ylabel('Degrees Latitude')
        set(gca,'FontSize',18)
%         printFig(gcf, 'MACA_Percent_Difference_MAX', [14 14], 'png', 150)
    case 'mean'
        var_plot_mean = 100*(((nan_mean(spd_avg(:,:,years1),3)-nan_mean(spd_avg(:,:,years2),3))./nan_mean(spd_avg(:,:,years1),3)));
        pcolor(M.lon,M.lat,var_plot_mean)
        shading flat
        c = colorbar;
        c.Label.String = 'Percent Difference - Average Wind Speed';
        caxis([-3 3]);
        title('Difference in Average Wind Speeds: 1950 - 2020 vs 2020 - 2100')
        hold on
        plot(wa_lon, wa_lat,'Color','k')
        axis equal
        xlim([min(M.lon(:)) -122.1])
        xlabel('Degrees Longitude')
        ylabel('Degrees Latitude')
        set(gca,'FontSize',18)
        printFig(gcf, 'MACA_Percent_Difference_MEAN', [14 14], 'png', 150)
end
%% Multi window percent difference plot 
clf
p_left = .03;
p_right = .05;
p_top = .08;
p_bot = .035;
p_spacing = .02;
p_wid = (1-p_right-p_left-p_spacing)/4;
p_height = (1-p_top-p_bot-p_spacing)/3;

metric = 'mean';

left_bound = min(M.lon(:));
right_bound = -122.01;
row = 2; 
col = 0;
mycolors = flipud(cbrewer('div','RdBu',N_c));
colormap(mycolors)
switch metric 
    case 'max'
        for x = 2:length(starts)
            axes('position',[p_left+col*(p_wid+p_spacing) p_bot+row*(p_height+p_spacing) p_wid p_height])
            var_plot_max = 100*((nan_mean(spd_max(:,:,dec_yrs{1,x-1}),3))-nan_mean(spd_max(:,:,dec_yrs{1,x}),3))./nan_mean(spd_max(:,:,dec_yrs{1,x-1}),3);
%             disp([min(var_plot_max(:)) max(var_plot_max(:))]);
            pcolor(M.lon,M.lat,var_plot_max) % plot
            set(gca,'FontSize', 8) % font size
            tit = sprintf('%d - %d vs %d - %d',starts(x-1),ends(x-1),starts(x),ends(x));
            title(tit)
            caxis([-3 3])
            xlim([left_bound right_bound])
            shading flat
            col = col + 1;
            if x == 2 || x == 6
                set(gca,'XTickLabel',[])
            elseif x == 11 || x == 12 || x == 13
                set(gca,'YTickLabel',[])
            elseif x == 3 || x == 4 || x == 5 || x == 6 || x == 7 || x == 8 || x == 9
                set(gca,'XTickLabel',[])
                set(gca,'YTickLabel',[])
            end
            if col > 3
                col = 0;
                row = row - 1;
            end
            hold on
            plot(wa_lon, wa_lat,'Color','k')%[.7 .7 .7]) % Add WA coastline
            pause
        end
        chan = colorbar('Location','NorthOutSide','FontSize',8);
        set(chan,'Position',[p_left 1-.04 ((p_left + 4*p_wid + 3*p_spacing)-p_left) .02])
        printFig(gcf, 'MACA_Percent_Difference_MAX_MULTIWINDOW', [14 14], 'png', 150)
    case 'mean'
        for x = 2:length(starts)
            axes('position',[p_left+col*(p_wid+p_spacing) p_bot+row*(p_height+p_spacing) p_wid p_height])
            var_plot_max = 100*((nan_mean(spd_avg(:,:,dec_yrs{1,x-1}),3))-nan_mean(spd_avg(:,:,dec_yrs{1,x}),3))./nan_mean(spd_avg(:,:,dec_yrs{1,x-1}),3);
%             var_plot_max = 100*(((nan_mean(M.spd(:,:,i{1,x}),3)-nan_mean(M.spd(:,:,i{x-1}),3))./nan_mean(M.spd(:,:,i{1,x}),3)));
            disp([min(var_plot_max(:)) max(var_plot_max(:))]);
            pcolor(M.lon,M.lat,var_plot_max) % plot
            set(gca,'FontSize', 8) % font size
            tit = sprintf('%d - %d vs %d - %d',starts(x-1),ends(x-1),starts(x),ends(x));
            title(tit)
            caxis([-1.5 1.5])
            xlim([left_bound right_bound])
            shading flat
            col = col + 1;
            if x == 2 || x == 6
                set(gca,'XTickLabel',[])
            elseif x == 11 || x == 12 || x == 13
                set(gca,'YTickLabel',[])
            elseif x == 3 || x == 4 || x == 5 || x == 6 || x == 7 || x == 8 || x == 9
                set(gca,'XTickLabel',[])
                set(gca,'YTickLabel',[])
            end
            if col > 3
                col = 0;
                row = row - 1;
            end
            hold on
            plot(wa_lon, wa_lat,'Color','k')%,[.7 .7 .7]) % Add WA coastline
            pause
        end
        chan = colorbar('Location','NorthOutSide','FontSize',10);
        set(chan,'Position',[p_left 1-.04 ((p_left + 4*p_wid + 3*p_spacing)-p_left) .02])
        printFig(gcf, 'MACA_Percent_Difference_MEAN_MULTIWINDOW', [14 14], 'png', 150)
end
%% Find Averages through time - Line Plot

% % Find Indices greater than determined speed
% spd_val = 15;
% t1_spd = M.spd(:,:,period1); t1_spd = t1_spd(:); % Grab all of the times I want and vectorize speed and direction
% t2_spd = M.spd(:,:,period2); t2_spd = t2_spd(:);
% t1_dir = M.dir(:,:,period1); t1_dir = t1_dir(:);
% t2_dir = M.dir(:,:,period2); t2_dir = t2_dir(:);
% t1_i = find(t1_spd >= spd_val); % find indices above specified threshold 
% t2_i = find(t2_spd >= spd_val); 

clf
% Plot the averages
a = 1:13;
b = zeros(1,13);
for x = 1:13
    temp_spds = spd_avg(:,:,dec_yrs{1,x}); temp_spds = temp_spds(:);
    b(x) = nan_mean(temp_spds);
end
plot(a,b)
hold on 
plot(a,b,'bo','MarkerFaceColor','b')
% use function and calculate R-squared
p = polyfit(a,b,1); % P1 is slope, P2 is intercept
f = polyval(p,a); 
[r2, rmse] = rsquare(b,f); 
hold on
l = plot(a,f,'--','Color',[.7 .7 .7],'LineWidth',.7);
% Play with xaxis labels
set(gca,'xtick',1:13)
xl = {'1950-1980','1960-1990','1970-2000','1980-2010','1990-2020',...
    '2000-2030','2010-2040','2020-2050','2030-2060','2040-2070','2050-2080','2060-2090','2070-2100'};
set(gca,'xticklabel',xl)
set(gca,'fontsize',8)
xlabel('Time - 30 year segments')
ylabel('Average Windspeed [m/s]')
grid on
% add slope to figure 
txt1 = sprintf('Slope: %4.2f',p(1));
text(1, 3.7, txt1,'FontSize',20)
printFig(gcf, 'Trend_WNDSPD_MEAN_MACA', [14 14], 'png', 150)

%% Find Maxes through time - Line Plot
clf
% Plot the averages
a = 1:13;
b = zeros(1,13);
for x = 1:13
    temp_spds = spd_max(:,:,dec_yrs{1,x}); temp_spds = temp_spds(:);
    b(x) = nan_mean(temp_spds);
end
plot(a,b)
hold on 
plot(a,b,'bo','MarkerFaceColor','b')
% use function and calculate R-squared
p = polyfit(a,b,1); % P1 is slope, P2 is intercept
f = polyval(p,a); 
[r2, rmse] = rsquare(b,f); 
hold on
l = plot(a,f,'--','Color',[.7 .7 .7],'LineWidth',.7);
% Play with xaxis labels
set(gca,'xtick',1:13)
xl = {'1950-1980','1960-1990','1970-2000','1980-2010','1990-2020',...
    '2000-2030','2010-2040','2020-2050','2030-2060','2040-2070','2050-2080','2060-2090','2070-2100'};
set(gca,'xticklabel',xl)
set(gca,'fontsize',8)
xlabel('Time - 30 year segments')
ylabel('Average Maximum Windspeed [m/s]')
grid on
% add slope to figure 
txt1 = sprintf('Slope: %4.3f',p(1));
text(1, 10.3, txt1,'FontSize',20)
printFig(gcf, 'Trend_WNDSPD_MAX_MACA', [14 14], 'png', 150)

%% Bin Winds by Magnitude
% Find Indices greater than determined speed
spd_val = 20;
t1_spd = M.spd(:,:,period1); t1_spd = t1_spd(:); % Grab all of the times I want and vectorize speed and direction
t2_spd = M.spd(:,:,period2); t2_spd = t2_spd(:);
t1_i = find(t1_spd >= spd_val); % find indices above specified threshold 
t2_i = find(t2_spd >= spd_val); 

clf
x_locs = 20:.2:23;
x_vec = 1:1:length(x_locs)-1;
h1 = histogram(t1_spd(t1_i),x_locs,'Normalization','probability');
l1 = h1.Values;
hold on
h2 = histogram(t2_spd(t2_i),x_locs,'Normalization','probability','FaceAlpha',0.4);
l2 = h2.Values;
legend('1950 - 2020','2020 - 2100')
xlabel('Wind Speed [m/s]')
ylabel('Probability')
tit = sprintf('MACA_Distribtuions_Binned_Wind_%d',spd_val);
printFig(gcf, tit, [14 14], 'png', 150)

%% Bin Winds by Direction
nor_i1 = find(M.dir(:,:,period1) >= 315 | M.dir(:,:,period1) <= 45);
nor_i2 = find(M.dir(:,:,period2) >= 315 | M.dir(:,:,period2) <= 45);
sou_i1 = find(M.dir(:,:,period1) >= 135 & M.dir(:,:,period1) <= 225);
sou_i2 = find(M.dir(:,:,period2) >= 135 & M.dir(:,:,period2) <= 225);
east_i1 = find(M.dir(:,:,period1) >= 45 & M.dir(:,:,period1) <=  135);
east_i2 = find(M.dir(:,:,period2) >= 45 & M.dir(:,:,period2) <=  135);
west_i1 = find(M.dir(:,:,period1) >= 225 & M.dir(:,:,period1) <=  315);
west_i2 = find(M.dir(:,:,period2) >= 225 & M.dir(:,:,period2) <=  315);

% Distribtuions for North Winds and changing through time
clf
x_locs = 1:1:16;
x_vec = 1:1:length(x_locs)-1;
h1 = histogram(t1_spd(nor_i1),x_locs,'Normalization','probability');
l1 = h1.Values;
hold on
h2 = histogram(t2_spd(nor_i2),x_locs,'Normalization','probability','FaceAlpha',0.4);
l2 = h2.Values;
legend('1950 - 2020','2020 - 2100')
xlabel('Wind Speed [m/s]')
ylabel('Probability')
title('North Winds [315\circ - 45\circ]')
tit = sprintf('MACA_Distribtuions_North_Binned_Wind');
printFig(gcf, tit, [14 14], 'png', 150)


% Distribtuions for South winds changing through time
clf
x_locs = 1:1:15;
x_vec = 1:1:length(x_locs)-1;
h1 = histogram(t1_spd(sou_i1),x_locs,'Normalization','probability');
l1 = h1.Values;
hold on
h2 = histogram(t2_spd(sou_i2),x_locs,'Normalization','probability','FaceAlpha',0.4);
l2 = h2.Values;
legend('1950 - 2020','2020 - 2100')
xlabel('Wind Speed [m/s]')
ylabel('Probability')
title('South Winds [135\circ - 225\circ]')
tit = sprintf('MACA_Distribtuions_South_Binned_Wind');
printFig(gcf, tit, [14 14], 'png', 150)

% Distribtuions for East winds changing through time
clf
x_locs = 1:1:15;
x_vec = 1:1:length(x_locs)-1;
h1 = histogram(t1_spd(east_i1),x_locs,'Normalization','probability');
l1 = h1.Values;
hold on
h2 = histogram(t2_spd(east_i2),x_locs,'Normalization','probability','FaceAlpha',0.4);
l2 = h2.Values;
legend('1950 - 2020','2020 - 2100')
xlabel('Wind Speed [m/s]')
ylabel('Probability')
title('East Winds [45\circ - 135\circ]')
tit = sprintf('MACA_Distribtuions_East_Binned_Wind');
printFig(gcf, tit, [14 14], 'png', 150)

% Distribtuions for West winds changing through time
clf
x_locs = 1:1:15;
x_vec = 1:1:length(x_locs)-1;
h1 = histogram(t1_spd(west_i1),x_locs,'Normalization','probability');
l1 = h1.Values;
hold on
h2 = histogram(t2_spd(west_i2),x_locs,'Normalization','probability','FaceAlpha',0.4);
l2 = h2.Values;
legend('1950 - 2020','2020 - 2100')
xlabel('Wind Speed [m/s]')
ylabel('Probability')
title('West Winds [225\circ - 315\circ]')
tit = sprintf('MACA_Distribtuions_West_Binned_Wind');
printFig(gcf, tit, [14 14], 'png', 150)
%% Heat Map of direction 
colormap(hsv)

% Time period 1
clf
pcolor(M.lon, M.lat, nan_mean(dir_avg(:,:,years1),3))
shading flat
xlabel('Degrees Longitude')
ylabel('Degrees Latitude')
c = colorbar;
c.Label.String = 'Degrees';
title('Average Wind Direction [degrees]: 1950 - 2020')
caxis([0 360])
hold on 
plot(wa_lon, wa_lat,'k')
axis equal
xlim([min(M.lon(:)) -122.1])
set(gca,'FontSize',18)
printFig(gcf, 'MACA_Average_DIR_1950_2020', [8.5 11], 'png', 150)

% Time period 2
clf
pcolor(M.lon, M.lat, nan_mean(dir_avg(:,:,years2),3))
shading flat
xlabel('Degrees Longitude')
ylabel('Degrees Latitude')
c = colorbar;
c.Label.String = 'Degrees';
title('Wind Direction [degrees]: 2020 - 2100')
caxis([0 360])
hold on 
plot(wa_lon, wa_lat,'k')
axis equal
xlim([min(M.lon(:)) -122.1])
set(gca,'FontSize',18)
printFig(gcf, 'MACA_Average_DIR_2020_2100', [8.5 11], 'png', 150)

% % Percent difference between the two
% clf
% N_c = 100;
% mycolors = flipud(cbrewer('div','RdBu',N_c));
% colormap(mycolors)
% var_plot_dir = 100*(((nan_mean(dir_avg(:,:,years1),3))-nan_mean(dir_avg(:,:,years2),3))./nan_mean(dir_avg(:,:,years1),3));
% pcolor(M.lon,M.lat,var_plot_dir)
% shading flat
% c = colorbar;
% c.Label.String = 'Percent Difference - Average Wind Direction ';
% caxis([-2 2]);
% title('Difference in Average Wind Direction : 1950 - 2020 vs 2020 - 2100')
% hold on
% plot(wa_lon, wa_lat,'Color','k')
% axis equal
% xlim([min(M.lon(:)) -122.1])
% xlabel('Degrees Longitude')
% ylabel('Degrees Latitude')
% printFig(gcf, 'MACA_Percent_Difference_Direction_AVG', [14 14], 'png', 150)

% Degree difference between the two
clf
N_c = 100;
mycolors = flipud(cbrewer('div','RdBu',N_c));
colormap(mycolors)
var_plot_deg = (nan_mean(dir_avg(:,:,years1),3))-nan_mean(dir_avg(:,:,years2),3);
pcolor(M.lon,M.lat,var_plot_deg)
shading flat
c = colorbar;
c.Label.String = 'Difference in Direction [degree]';
caxis([-15 15]);
title('Difference in Average Wind Direction : 1950 - 2020 vs 2020 - 2100')
hold on
plot(wa_lon, wa_lat,'Color','k')
axis equal
xlim([min(M.lon(:)) -122.1])
xlabel('Degrees Longitude')
ylabel('Degrees Latitude')
set(gca,'FontSize',18)
printFig(gcf, 'MACA_Difference_Degree_AVG', [14 14], 'png', 150)

% Percent difference between Max windspeed directions 
% clf
% N_c = 100;
% mycolors = flipud(cbrewer('div','RdBu',N_c));
% colormap(mycolors)
% var_plot_dir = 100*(((nan_mean(dir_max_mean(:,:,years1),3))-nan_mean(dir_max_mean(:,:,years2),3))./nan_mean(dir_max_mean(:,:,years1),3));
% pcolor(M.lon,M.lat,var_plot_dir)
% shading flat
% c = colorbar;
% c.Label.String = 'Percent Difference - Average Wind Direction ';
% caxis([-6 6]);
% title('Difference in Max Wind Direction : 1950 - 2020 vs 2020 - 2100')
% hold on
% plot(wa_lon, wa_lat,'Color','k')
% axis equal
% xlim([min(M.lon(:)) -122.1])
% xlabel('Degrees Longitude')
% ylabel('Degrees Latitude')
% printFig(gcf, 'MACA_Percent_Difference_Direction_MAX_SPEED_DIRECTION', [14 14], 'png', 150)

% Difference in degree of max wind speed average
clf
N_c = 100;
mycolors = flipud(cbrewer('div','RdBu',N_c));
colormap(mycolors)
var_plot_deg = ((nan_mean(dir_max_mean(:,:,years1),3))-nan_mean(dir_max_mean(:,:,years2),3));
pcolor(M.lon,M.lat,var_plot_deg)
shading flat
c = colorbar;
c.Label.String = 'Difference in Direction [degree]';
caxis([-30 30]);
title('Degree Difference in Max Wind Direction : 1950 - 2020 vs 2020 - 2100')
hold on
plot(wa_lon, wa_lat,'Color','k')
axis equal
xlim([min(M.lon(:)) -122.1])
xlabel('Degrees Longitude')
ylabel('Degrees Latitude')
set(gca,'FontSize',18)
printFig(gcf, 'MACA_Degree_Difference_MAX_SPEED_DIRECTION', [14 14], 'png', 150)



%% Multi Window Percent Difference Plot for Direction - Ignore for now
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % clf
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % p_left = .03;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % p_right = .05;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % p_top = .08;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % p_bot = .035;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % p_spacing = .02;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % p_wid = (1-p_right-p_left-p_spacing)/4;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % p_height = (1-p_top-p_bot-p_spacing)/3;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % left_bound = min(M.lon(:));
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % right_bound = -122.01;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % row = 2; 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % col = 0;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % N_c = 100;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % mycolors = flipud(cbrewer('div','RdBu',N_c));
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % colormap(mycolors)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % for x = 2:length(starts)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     axes('position',[p_left+col*(p_wid+p_spacing) p_bot+row*(p_height+p_spacing) p_wid p_height])
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     var_plot_dir = 100*(((nan_mean(M.dir(:,:,i{1,x}),3)-nan_mean(M.dir(:,:,i{x-1}),3))./nan_mean(M.dir(:,:,i{1,x}),3)));
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     pcolor(M.lon,M.lat,var_plot_dir) % plot
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     set(gca,'FontSize', 8) % font size
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     tit = sprintf('%d - %d vs %d - %d',starts(x-1),ends(x-1),starts(x),ends(x));
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     title(tit)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     caxis([-2 2])
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     xlim([left_bound right_bound])
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     shading flat
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     col = col + 1;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     if x == 2 || x == 6
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         set(gca,'XTickLabel',[])
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     elseif x == 11 || x == 12 || x == 13
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         set(gca,'YTickLabel',[])
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     elseif x == 3 || x == 4 || x == 5 || x == 6 || x == 7 || x == 8 || x == 9
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         set(gca,'XTickLabel',[])
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         set(gca,'YTickLabel',[])
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     if col > 3
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         col = 0;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         row = row - 1;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     hold on
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     plot(wa_lon, wa_lat,'Color','k')%[.7 .7 .7]) % Add WA coastline
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % chan = colorbar('Location','NorthOutSide','FontSize',8);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % set(chan,'Position',[p_left 1-.04 ((p_left + 4*p_wid + 3*p_spacing)-p_left) .02])
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % printFig(gcf, 'MACA_Percent_Difference_DIR_MULTIWINDOW', [14 14], 'png', 150)
%% Multi Window Difference Plot for Future Compared to Current
clf
p_left = .03;
p_right = .05;
p_top = .08;
p_bot = .035;
p_spacing = .02;
p_wid = (1-p_right-p_left-p_spacing)/4;
p_height = (1-p_top-p_bot-p_spacing)/2;

left_bound = min(M.lon(:));
right_bound = -122.01;

row = 1; 
col = 0;
N_c = 100;
mycolors = flipud(cbrewer('div','RdBu',N_c));
colormap(mycolors)
clear i
starts = 2030:10:2090;
ends = 2040:10:2100;
for xx = 1:length(starts)
    i{xx} = find(years >= starts(xx) & years <= ends(xx));
end


% metric = 'maxspd';
% metric = 'maxdir';
% metric = 'avgspd';
 metric = 'avgdir';

switch metric 
    case 'maxspd'
        for x = 1:length(starts)
            if x <= 4
                axes('position',[p_left+col*(p_wid+p_spacing) p_bot+row*(p_height+p_spacing) p_wid p_height])
            else
                axes('position',[p_left+(.5*p_wid)+col*(p_wid+p_spacing) p_bot+row*(p_height+p_spacing) p_wid p_height])
            end
            var_plot_spd = 100*(((nan_mean(spd_max(:,:,years1),3)-nan_mean(spd_max(:,:,i{1,x}),3))./nan_mean(spd_max(:,:,years1),3)));
            disp([min(var_plot_spd(:)) max(var_plot_spd(:))])
            pcolor(M.lon,M.lat,var_plot_spd) % plot
            shading flat
            set(gca,'FontSize', 8) % font size
            tit = sprintf('1950 - 2020 Max Speed vs %d Max Speed',starts(x));
            title(tit)
            caxis([-5 5])
            xlim([left_bound right_bound])
            shading flat
            col = col + 1;
            if x == 1
                set(gca,'XTickLabel',[])
            elseif x == 6 || x == 7 || x == 8
                set(gca,'YTickLabel',[])
            elseif x == 2 || x == 3 || x == 4
                set(gca,'XTickLabel',[])
                set(gca,'YTickLabel',[])
            end
            if col > 3
                col = 0;
                row = row - 1;
            end
            hold on
            plot(wa_lon, wa_lat,'Color','k')%[.7 .7 .7]) % Add WA coastline
            
        end
        chan = colorbar('Location','NorthOutSide','FontSize',8);
        chan.Label.String = 'Percent Difference';
        set(chan,'Position',[p_left 1-.059 ((p_left + 4*p_wid + 3*p_spacing)-p_left) .02])
        printFig(gcf, 'MACA_Percent_Difference_Speed_CurMAX2Future', [11 8.5], 'png', 150)
    case 'avgspd'
        for x = 1:length(starts)
            if x <= 4
                axes('position',[p_left+col*(p_wid+p_spacing) p_bot+row*(p_height+p_spacing) p_wid p_height])
            else
                axes('position',[p_left+(.5*p_wid)+col*(p_wid+p_spacing) p_bot+row*(p_height+p_spacing) p_wid p_height])
            end
            var_plot_spd = 100*(((nan_mean(spd_avg(:,:,years1),3)-nan_mean(spd_avg(:,:,i{1,x}),3))./nan_mean(spd_avg(:,:,years1),3)));
            %disp([min(var_plot_spd(:)) max(var_plot_spd(:))])
            pcolor(M.lon,M.lat,var_plot_spd) % plot
            shading flat
            set(gca,'FontSize', 8) % font size
            tit = sprintf('1950 - 2020 Avg Speed vs %d Avg Speed',starts(x));
            title(tit)
            caxis([-4 4])
            xlim([left_bound right_bound])
            shading flat
            col = col + 1;
            if x == 1
                set(gca,'XTickLabel',[])
            elseif x == 6 || x == 7 || x == 8
                set(gca,'YTickLabel',[])
            elseif x == 2 || x == 3 || x == 4
                set(gca,'XTickLabel',[])
                set(gca,'YTickLabel',[])
            end
            if col > 3
                col = 0;
                row = row - 1;
            end
            hold on
            plot(wa_lon, wa_lat,'Color','k')%[.7 .7 .7]) % Add WA coastline
            
        end
        chan = colorbar('Location','NorthOutSide','FontSize',8);
        chan.Label.String = 'Percent Difference';
        set(chan,'Position',[p_left 1-.059 ((p_left + 4*p_wid + 3*p_spacing)-p_left) .02])
        printFig(gcf, 'MACA_Percent_Difference_Speed_CurAvg2Future', [11 8.5], 'png', 150)
    case 'maxdir'
        for x = 1:length(starts)
            if x <= 4
                axes('position',[p_left+col*(p_wid+p_spacing) p_bot+row*(p_height+p_spacing) p_wid p_height])
            else
                axes('position',[p_left+(.5*p_wid)+col*(p_wid+p_spacing) p_bot+row*(p_height+p_spacing) p_wid p_height])
            end
            var_plot_dir = (nan_mean(dir_max_mean(:,:,years1),3)-nan_mean(dir_max_mean(:,:,i{1,x}),3));
            disp([min(var_plot_dir(:)) max(var_plot_dir(:))])
            pcolor(M.lon,M.lat,var_plot_dir) % plot
            shading flat
            set(gca,'FontSize', 8) % font size
            tit = sprintf('1950 - 2020 "Max" Direction vs %d "Max" Direction',starts(x));
            title(tit)
            caxis([-60 60])
            xlim([left_bound right_bound])
            shading flat
            col = col + 1;
            if x == 1
                set(gca,'XTickLabel',[])
            elseif x == 6 || x == 7 || x == 8
                set(gca,'YTickLabel',[])
            elseif x == 2 || x == 3 || x == 4
                set(gca,'XTickLabel',[])
                set(gca,'YTickLabel',[])
            end
            if col > 3
                col = 0;
                row = row - 1;
            end
            hold on
            plot(wa_lon, wa_lat,'Color','k')%[.7 .7 .7]) % Add WA coastline
            
        end
        chan = colorbar('Location','NorthOutSide','FontSize',8);
        chan.Label.String = 'Degree Difference';
        set(chan,'Position',[p_left 1-.059 ((p_left + 4*p_wid + 3*p_spacing)-p_left) .02])
        printFig(gcf, 'MACA_Degree_Difference_Direction_CurMAX2Future', [11 8.5], 'png', 150)
    case 'avgdir'
        for x = 1:length(starts)
            if x <= 4
                axes('position',[p_left+col*(p_wid+p_spacing) p_bot+row*(p_height+p_spacing) p_wid p_height])
            else
                axes('position',[p_left+(.5*p_wid)+col*(p_wid+p_spacing) p_bot+row*(p_height+p_spacing) p_wid p_height])
            end
            var_plot_dir = (nan_mean(dir_avg(:,:,years1),3)-nan_mean(dir_avg(:,:,i{1,x}),3));
            disp([min(var_plot_dir(:)) max(var_plot_dir(:))])
            pcolor(M.lon,M.lat,var_plot_dir) % plot
            shading flat
            set(gca,'FontSize', 8) % font size
            tit = sprintf('1950 - 2020 Avg Direction vs %d Avg Direction',starts(x));
            title(tit)
            caxis([-30 30])
            xlim([left_bound right_bound])
            shading flat
            col = col + 1;
            if x == 1
                set(gca,'XTickLabel',[])
            elseif x == 6 || x == 7 || x == 8
                set(gca,'YTickLabel',[])
            elseif x == 2 || x == 3 || x == 4
                set(gca,'XTickLabel',[])
                set(gca,'YTickLabel',[])
            end
            if col > 3
                col = 0;
                row = row - 1;
            end
            hold on
            plot(wa_lon, wa_lat,'Color','k')%[.7 .7 .7]) % Add WA coastline
            
        end
        chan = colorbar('Location','NorthOutSide','FontSize',8);
        chan.Label.String = 'Degree Difference';
        set(chan,'Position',[p_left 1-.059 ((p_left + 4*p_wid + 3*p_spacing)-p_left) .02])
        printFig(gcf, 'MACA_Degree_Difference_Direction_CurAVG2Future', [11 8.5], 'png', 150)
end
%% 90th Percentile of winds 
clf
p_left = .03;
p_right = .05;
p_top = .08;
p_bot = .035;
p_spacing = .02;
p_wid = (1-p_right-p_left-p_spacing)/4;
p_height = (1-p_top-p_bot-p_spacing)/2;

left_bound = min(M.lon(:));
right_bound = -122.01;

row = 1; 
col = 0;
N_c = 100;
mycolors = flipud(cbrewer('div','RdBu',N_c));
colormap(mycolors)
clear i
starts = 2030:10:2090;
ends = 2040:10:2100;
for xx = 1:length(starts)
    i{xx} = find(years >= starts(xx) & years <= ends(xx));
end
pct = 99;
for x = 1:length(starts)
    if x <= 4
        axes('position',[p_left+col*(p_wid+p_spacing) p_bot+row*(p_height+p_spacing) p_wid p_height])
    else
        axes('position',[p_left+(.5*p_wid)+col*(p_wid+p_spacing) p_bot+row*(p_height+p_spacing) p_wid p_height])
    end
    var_plot_pct = 100*(((prctile(M.spd(:,:,period1),pct,3)-prctile(M.spd(:,:,i{1,x}),pct,3))./prctile(M.spd(:,:,period1),pct,3)));
    %disp([min(var_plot_pct(:)) max(var_plot_pct(:))])
    pcolor(M.lon,M.lat,var_plot_pct) % plot
    shading flat
    set(gca,'FontSize', 8) % font size
    tit = sprintf('1950 - 2020: %d^{th} Pct vs %d %d^{th} Pct',pct,starts(x),pct);
    title(tit)
    caxis([-80 80]) % -50 and 50 for 90, -80 and 80 for 99th
    xlim([left_bound right_bound])
    shading flat
    col = col + 1;
    if x == 1 
        set(gca,'XTickLabel',[])
    elseif x == 6 || x == 7 || x == 8
        set(gca,'YTickLabel',[])
    elseif x == 2 || x == 3 || x == 4 
        set(gca,'XTickLabel',[])
        set(gca,'YTickLabel',[])
    end
    if col > 3
        col = 0;
        row = row - 1;
    end
    hold on
    plot(wa_lon, wa_lat,'Color','k')%[.7 .7 .7]) % Add WA coastline
    
end
chan = colorbar('Location','NorthOutSide','FontSize',8);
chan.Label.String = 'Percent Difference';
set(chan,'Position',[p_left 1-.059 ((p_left + 4*p_wid + 3*p_spacing)-p_left) .02])
tit = sprintf('MACA_Percent_Difference_%d_Pctile_SPD',pct);
printFig(gcf, tit, [11 8.5], 'png', 150)

%% QQ Plot for MACA 
hist = find(year(M.time) >= 1950 & year(M.time) <= 2020);
fut = find(year(M.time) >= 2030 & year(M.time) <= 2100);
% Grab indices and get equal times for QQ plots 
hist(1:50) = [];

metric = 'spd';
switch metric 
    case 'spd'
        data_hist = M.spd(:,:,hist); data_hist = data_hist(:);
        data_fut = M.spd(:,:,fut); data_fut = data_fut(:);
        fname = 'MACA_85_QQ_Speed';
    case 'dir'
        data_hist = M.dir(:,:,hist);
        data_fut = M.dir(:,:,fut);
        fname = 'MACA_85_QQ_Direction'
end
cdf_hist = sort(data_hist,'ascend');
cdf_fut = sort(data_fut,'ascend');
x_axis = linspace(0,1,length(cdf_hist));
clf
plot(cdf_hist,cdf_fut,'LineWidth',2)
grid on
hold on 
line([0 50],[0 50],'Color','k','LineStyle','--')
xlabel('MACA Historic: 1950 - 2020');
ylabel('MACA Future: 2030 - 2100');
xlim([0 25])
ylim([0 25])

set(gca,'FontSize',14)

printFig(gcf,'MACA_85_QQ_Speed',[14 14],'png',300)



