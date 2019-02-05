% Loads in MACA data and creates Heat Maps for Change Analysis

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
scenario = 'rcp45';
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
        
        Mh_u = ncread('E:\Abbas\Model_Met_Forcings\MACA\maca_u_GFDL_ESM2M_historical_1950_2005.nc','eastward_wind'); % THIS IS FOR PLOTTING LATER
        
end

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

max2grab = 1; 
% number of maxes to grab for average max per year
% Grab n-max values per year for each grid cell, and then average them to
% represent the 'max' per year per grid cell.  Also calculate the average
% direction for the maxes.  


%%%
% U AND V Shit 
%     wndspd(tt,:) = sqrt(u10.^2 + v10.^2);
%     wnddir_temp = (180/pi)*atan2(v10,u10);
%     
%     % Rotate winds to compass directions
%     wnddir_temp = 90 - wnddir_temp;
%     wnddir_temp(wnddir_temp<0)=wnddir_temp(wnddir_temp<0)+360;
%     
%%%


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
fname = sprintf('MACA_Stats_Data_%s_MaxsGrabbed_%d',scenario,max2grab);
lat = M.lat; lon = M.lon;
save(fname,'spd_max','spd_avg','years1','years2','years','dir_max_mean','dir_avg','u_max','v_max','u_avg','v_avg','spd_90','spd_95','spd_99','lat','lon');
% load('MACA_Stats_Data_rcp85_MaxsGrabbed_1.mat');
return
%% Plots Heat map of Longer period differences 

metric = 'mean';
param = 'absolute';

clf
N_c = 100;
mycolors = flipud(cbrewer('div','RdBu',N_c));
colormap(mycolors)

switch metric 
    case 'max'
        switch param
            case 'relative'
                var_plot_max = 100*(((nan_mean(spd_max(:,:,years2),3)-nan_mean(spd_max(:,:,years1),3))./nan_mean(spd_max(:,:,years2),3)));
                pcolor(lon,lat,var_plot_max)
                shading interp
                c = colorbar;
                c.Label.String = 'Percent Difference in Wind Speed [%]';
                caxis([-5 5]);
                fname = 'MACA_Percent_Difference_MAX';
            case 'absolute'
                var_plot_max = nan_mean(spd_max(:,:,years2),3)-nan_mean(spd_max(:,:,years1),3);
%                 var_plot_max = 100*(((nan_mean(spd_max(:,:,years1),3)-nan_mean(spd_max(:,:,years2),3))./nan_mean(spd_max(:,:,years1),3)));
                pcolor(lon,lat,var_plot_max)
                shading interp
                c = colorbar;
                c.Label.String = 'Difference in Wind Speed [m/s]';
                caxis([-.5 .5]);
                fname = sprintf('MACA_Absolute_Difference_MAX_%s',scenario);
        end
        title('Difference in Mean Yearly Max Wind Speeds: 1950 - 2020 vs 2020 - 2100')
        hold on
        plot(wa_lon, wa_lat,'Color','k')
        axis equal
        xlim([min(lon(:)) -122.1])
%         ylim([47 49.2])
        xlabel('Degrees Longitude')
        ylabel('Degrees Latitude')
        set(gca,'FontSize',18)
        printFig(gcf, fname, [14 14], 'png', 150)
    case 'mean'
        switch param
            case 'relative'
                var_plot_mean = 100*(((nan_mean(spd_avg(:,:,years1),3)-nan_mean(spd_avg(:,:,years2),3))./nan_mean(spd_avg(:,:,years1),3)));
                pcolor(lon,lat,var_plot_mean)
                shading interp
                c = colorbar;
                c.Label.String = 'Percent Difference - Average Wind Speed';
                caxis([-3 3]);
                fname = 'MACA_Percent_Difference_AVERAGE';
            case 'absolute'
                var_plot_mean = nan_mean(spd_avg(:,:,years2),3)-nan_mean(spd_avg(:,:,years1),3);
%                 var_plot_max = 100*(((nan_mean(spd_max(:,:,years1),3)-nan_mean(spd_max(:,:,years2),3))./nan_mean(spd_max(:,:,years1),3)));
                pcolor(lon,lat,var_plot_mean)
                shading interp
                c = colorbar;
                c.Label.String = 'Difference in Avg. Wind Speed [m/s]';
                caxis([-.2 .2]);
                fname = sprintf('MACA_Absolute_Difference_AVERAGE_%s',scenario);
        end
        title('Difference in Average Wind Speeds: 1950 - 2020 vs 2020 - 2100')
        hold on
        plot(wa_lon, wa_lat,'Color','k')
        axis equal
        xlim([min(lon(:)) -122.1])
        xlabel('Degrees Longitude')
        ylabel('Degrees Latitude')
        set(gca,'FontSize',18)
        printFig(gcf, fname, [14 14], 'png', 150)
    case '90'
        switch param
            case 'relative'
                var_plot_90 = 100*(((nan_mean(spd_90(:,:,years1),3)-nan_mean(spd_90(:,:,years2),3))./nan_mean(spd_90(:,:,years1),3)));
                pcolor(lon,lat,var_plot_90)
                shading interp
                c = colorbar;
                c.Label.String = 'Percent Difference - Average Wind Speed';
                caxis([-3 3]);
                fname =  'MACA_Percent_Difference_90'
            case 'absolute'
                var_plot_90 = nan_mean(spd_90(:,:,years2),3)-nan_mean(spd_90(:,:,years1),3);
%                 var_plot_max = 100*(((nan_mean(spd_max(:,:,years1),3)-nan_mean(spd_max(:,:,years2),3))./nan_mean(spd_max(:,:,years1),3)));
                pcolor(lon,lat,var_plot_90)
                shading interp
                c = colorbar;
                c.Label.String = 'Difference in 90th Wind Speed [m/s]';
                caxis([-.31 .31]);
                fname =  'MACA_Absolute_Difference_90'
        end
        title('Difference in 90th Pct Wind Speeds: 1950 - 2020 vs 2020 - 2100')
        hold on
        plot(wa_lon, wa_lat,'Color','k')
        axis equal
        xlim([min(lon(:)) -122.1])
        xlabel('Degrees Longitude')
        ylabel('Degrees Latitude')
        set(gca,'FontSize',18)
        printFig(gcf,fname, [14 14], 'png', 150)
    case '95'
        switch param
            case 'relative'
                var_plot_95 = 100*(((nan_mean(spd_95(:,:,years1),3)-nan_mean(spd_95(:,:,years2),3))./nan_mean(spd_95(:,:,years1),3)));
                pcolor(lon,lat,var_plot_95)
                shading interp
                c = colorbar;
                c.Label.String = 'Percent Difference - Average Wind Speed';
                caxis([-3 3]);
                fname = 'MACA_Percent_Difference_95'
            case 'absolute'
                var_plot_95 = nan_mean(spd_95(:,:,years2),3)-nan_mean(spd_95(:,:,years1),3);
%                 var_plot_max = 100*(((nan_mean(spd_max(:,:,years1),3)-nan_mean(spd_max(:,:,years2),3))./nan_mean(spd_max(:,:,years1),3)));
                pcolor(lon,lat,var_plot_95)
                shading interp
                c = colorbar;
                c.Label.String = 'Difference in 95th Wind Speed [m/s]';
                caxis([-.31 .31]);
                fname = 'MACA_Absolute_Difference_95'
        end
        title('Difference in 95th Pct Wind Speeds: 1950 - 2020 vs 2020 - 2100')
        hold on
        plot(wa_lon, wa_lat,'Color','k')
        axis equal
        xlim([min(lon(:)) -122.1])
        xlabel('Degrees Longitude')
        ylabel('Degrees Latitude')
        set(gca,'FontSize',18)
        printFig(gcf, fname, [14 14], 'png', 150)
    case '99'
        switch param
            case 'relative'
                var_plot_99 = 100*(((nan_mean(spd_99(:,:,years2),3)-nan_mean(spd_99(:,:,years1),3))./nan_mean(spd_99(:,:,years2),3)));
                pcolor(lon,lat,var_plot_99)
                shading interp
                c = colorbar;
                c.Label.String = 'Percent Difference - Average Wind Speed';
                caxis([-3 3]);
                fname = 'MACA_Percent_Difference_99';
            case 'absolute'
                var_plot_99 = nan_mean(spd_99(:,:,years2),3)-nan_mean(spd_99(:,:,years1),3);
%                 var_plot_max = 100*(((nan_mean(spd_max(:,:,years1),3)-nan_mean(spd_max(:,:,years2),3))./nan_mean(spd_max(:,:,years1),3)));
                pcolor(lon,lat,var_plot_99)
                shading interp
                c = colorbar;
                c.Label.String = 'Difference in 99th Wind Speed [m/s]';
                caxis([-.4 .4]);
                fname = 'MACA_Absolute_Difference_99';
        end
        title('Difference in 99th Pct Wind Speeds: 1950 - 2020 vs 2020 - 2100')
        hold on
        plot(wa_lon, wa_lat,'Color','k')
        axis equal
        xlim([min(lon(:)) -122.1])
        xlabel('Degrees Longitude')
        ylabel('Degrees Latitude')
        set(gca,'FontSize',18)
        printFig(gcf, fname, [14 14], 'png', 150)
    case 'dirMax'
        switch param
            case 'relative'
                var_plot_99 = 100*(((nan_mean(dir_max_mean(:,:,years2),3)-nan_mean(dir_max_mean(:,:,years1),3))./nan_mean(dir_max_mean(:,:,years2),3)));
                pcolor(lon,lat,var_plot_99)
                shading interp
                c = colorbar;
                c.Label.String = 'Percent Difference - Average Wind Speed';
                caxis([-3 3]);
                fname = 'MACA_Percent_Difference_dirMAX';
                
            %%%%%%%%%%%%%%%%%% CVI    
            case 'absolute'
                var_plot_dirMax = circ_mean(dir_max_mean(:,:,years2),[],3)-circ_mean(dir_max_mean(:,:,years1),[],3);
%                 var_plot_max = 100*(((nan_mean(spd_max(:,:,years1),3)-nan_mean(spd_max(:,:,years2),3))./nan_mean(spd_max(:,:,years1),3)));
                pcolor(lon,lat,var_plot_dirMax)
                shading interp
                c = colorbar;
                c.Label.String = 'Difference in Wind Direction [degrees]';
                caxis([-5 5]);
                fname = 'MACA_Absolute_Difference_dirMAX';
        end
        title('Difference in Wind Direction for Max Speeds: 1950 - 2020 vs 2020 - 2100')
        hold on
        plot(wa_lon, wa_lat,'Color','k')
        axis equal
        xlim([min(lon(:)) -122.1])
        xlabel('Degrees Longitude')
        ylabel('Degrees Latitude')
        set(gca,'FontSize',18)
        printFig(gcf, fname, [14 14], 'png', 150)
    case 'dirAvg'
        switch param
            case 'relative'
                var_plot_99 = 100*(((nan_mean(dir_max_mean(:,:,years2),3)-nan_mean(dir_max_mean(:,:,years1),3))./nan_mean(dir_max_mean(:,:,years2),3)));
                pcolor(lon,lat,var_plot_99)
                shading interp
                c = colorbar;
                c.Label.String = 'Percent Difference - Average Wind Speed';
                caxis([-3 3]);
                fname = 'MACA_Percent_Difference_dirAVG';
            case 'absolute'
                var_plot_dirAvg = nan_mean(dir_avg(:,:,years2),3)-nan_mean(dir_avg(:,:,years1),3);
%                 var_plot_max = 100*(((nan_mean(spd_max(:,:,years1),3)-nan_mean(spd_max(:,:,years2),3))./nan_mean(spd_max(:,:,years1),3)));
                pcolor(lon,lat,var_plot_dirAvg)
                shading interp
                c = colorbar;
                c.Label.String = 'Difference in Wind Direction [degrees]';
                caxis([-15 15]);
                fname = 'MACA_Absolute_Difference_dirAVG';
        end
        title('Difference in Average Wind Direction: 1950 - 2020 vs 2020 - 2100')
        hold on
        plot(wa_lon, wa_lat,'Color','k')
        axis equal
        xlim([min(lon(:)) -122.1])
        xlabel('Degrees Longitude')
        ylabel('Degrees Latitude')
        set(gca,'FontSize',18)
        printFig(gcf, fname, [14 14], 'png', 150)
    case 'dirUV'
        switch param
            case 'relative'
                var_plot_99 = 100*(((nan_mean(dir_max_mean(:,:,years2),3)-nan_mean(dir_max_mean(:,:,years1),3))./nan_mean(dir_max_mean(:,:,years2),3)));
                pcolor(lon,lat,var_plot_99)
                shading interp
                c = colorbar;
                c.Label.String = 'Percent Difference - Average Wind Speed';
                caxis([-3 3]);
                fname = 'MACA_Percent_Difference_dirAVG';
            case 'absolute'
                var_plot_dirAvg = nan_mean(dir_avg(:,:,years2),3)-nan_mean(dir_avg(:,:,years1),3);
%                 var_plot_max = 100*(((nan_mean(spd_max(:,:,years1),3)-nan_mean(spd_max(:,:,years2),3))./nan_mean(spd_max(:,:,years1),3)));
                pcolor(lon,lat,var_plot_dirAvg)
                shading interp
                c = colorbar;
                c.Label.String = 'Difference in Wind Direction [degrees]';
                caxis([-15 15]);
                fname = sprintf('MACA_Absolute_Difference_dirAVG_%s',scenario);
        end
        title('Difference in Average Wind Direction: 1950 - 2020 vs 2020 - 2100')
        hold on
        plot(wa_lon, wa_lat,'Color','k')
        axis equal
        xlim([min(lon(:)) -122.1])
        xlabel('Degrees Longitude')
        ylabel('Degrees Latitude')
        set(gca,'FontSize',18)
        printFig(gcf, fname, [14 14], 'png', 150)
        
end


%% Multi Window Difference Plot for Future Compared to Current
clf
p_left = .038;
p_right = .05;
p_top = .1;
p_bot = .06;
p_spacing = .02;
p_wid = (1-p_right-p_left-p_spacing)/4;
p_height = (1-p_top-p_bot-p_spacing)/2;

left_bound = min(lon(:));
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

metric = 'avgspd';
% metric = 'maxspd';
% metric = 'maxdir';
% metric = 'avgspd';
%  metric = '99spd';

switch metric 
    case 'maxspd'
        for x = 1:length(starts)
            if x <= 4
                axes('position',[p_left+col*(p_wid+p_spacing) p_bot+row*(p_height+p_spacing) p_wid p_height])
            else
                axes('position',[p_left+(.5*p_wid)+col*(p_wid+p_spacing) p_bot+row*(p_height+p_spacing) p_wid p_height])
            end
            var_plot_spd = nan_mean(spd_max(:,:,i{1,x}),3)-nan_mean(spd_max(:,:,years1),3);
            pcolor(lon,lat,var_plot_spd) % plot
            shading interp
            set(gca,'FontSize', 8) % font size
            tit = sprintf('1950 - 2020 Max Speed vs %d Max Speed',starts(x));
            title(tit)
            caxis([-1 1])
            xlim([left_bound right_bound])
            col = col + 1;
            if x == 1
                set(gca,'XTickLabel',[])
                ylabel('Degrees Latitude')
            elseif x == 6 || x == 7 || x == 8
                set(gca,'YTickLabel',[])
                xlabel('Degrees Longitude')
            elseif x == 2 || x == 3 || x == 4
                set(gca,'XTickLabel',[])
                set(gca,'YTickLabel',[])
            else
                ylabel('Degrees Latitude')
                xlabel('Derees Longitude')
            end
            if col > 3
                col = 0;
                row = row - 1;
            end
            hold on
            plot(wa_lon, wa_lat,'Color','k')%[.7 .7 .7]) % Add WA coastline
            
        end
        chan = colorbar('Location','NorthOutSide','FontSize',8);
        chan.Label.String = 'Difference [m/s]';
        chan.FontSize = 11;
        set(chan,'Position',[p_left .925 ((p_left + 4*p_wid + 3*p_spacing)-p_left) .02])
        fname = sprintf('MACA_Absolute_Difference_Speed_CurMAX2Future_%s',scenario);
        printFig(gcf, fname, [13 8.5], 'png', 150)
    case '99spd'
        for x = 1:length(starts)
            if x <= 4
                axes('position',[p_left+col*(p_wid+p_spacing) p_bot+row*(p_height+p_spacing) p_wid p_height])
            else
                axes('position',[p_left+(.5*p_wid)+col*(p_wid+p_spacing) p_bot+row*(p_height+p_spacing) p_wid p_height])
            end
            var_plot_spd99 = nan_mean(spd_99(:,:,i{1,x}),3)-nan_mean(spd_99(:,:,years1),3);
            pcolor(lon,lat,var_plot_spd99) % plot
            shading interp
            set(gca,'FontSize', 8) % font size
            tit = sprintf('1950 - 2020 Max Speed vs %d Max Speed',starts(x));
            title(tit)
            caxis([-1 1])
            xlim([left_bound right_bound])
            col = col + 1;
            if x == 1
                set(gca,'XTickLabel',[])
                ylabel('Degrees Latitude')
            elseif x == 6 || x == 7 || x == 8
                set(gca,'YTickLabel',[])
                xlabel('Degrees Longitude')
            elseif x == 2 || x == 3 || x == 4
                set(gca,'XTickLabel',[])
                set(gca,'YTickLabel',[])
            else
                ylabel('Degrees Latitude')
                xlabel('Derees Longitude')
            end
            if col > 3
                col = 0;
                row = row - 1;
            end
            hold on
            plot(wa_lon, wa_lat,'Color','k')%[.7 .7 .7]) % Add WA coastline
            
        end
        chan = colorbar('Location','NorthOutSide','FontSize',8);
        chan.Label.String = 'Difference [m/s]';
        chan.FontSize = 11;
        set(chan,'Position',[p_left .925 ((p_left + 4*p_wid + 3*p_spacing)-p_left) .02])
        printFig(gcf, 'MACA_Absolute_Difference_Speed99Pctile_CurMAX2Future', [13 8.5], 'png', 150)
    case 'avgspd'
        for x = 1:length(starts)
            if x <= 4
                axes('position',[p_left+col*(p_wid+p_spacing) p_bot+row*(p_height+p_spacing) p_wid p_height])
            else
                axes('position',[p_left+(.5*p_wid)+col*(p_wid+p_spacing) p_bot+row*(p_height+p_spacing) p_wid p_height])
            end
            var_plot_spd = nan_mean(spd_avg(:,:,i{1,x}),3) - nan_mean(spd_avg(:,:,years1),3);
            %disp([min(var_plot_spd(:)) max(var_plot_spd(:))])
            pcolor(lon,lat,var_plot_spd) % plot
            shading flat
            set(gca,'FontSize', 8) % font size
            tit = sprintf('1950 - 2020 Avg Speed vs %d Avg Speed',starts(x));
            title(tit)
            caxis([-.3 .3])
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
        chan.FontSize = 11;
        set(chan,'Position',[p_left .925 ((p_left + 4*p_wid + 3*p_spacing)-p_left) .02])
        fname = sprintf('MACA_Percent_Difference_Speed_CurAvg2Future_%s',scenario);
        printFig(gcf, 'MACA_Percent_Difference_Speed_CurAvg2Future', [11 8.5], 'png', 150)
    case 'maxdir'
        for x = 1:length(starts)
            if x <= 4
                axes('position',[p_left+col*(p_wid+p_spacing) p_bot+row*(p_height+p_spacing) p_wid p_height])
            else
                axes('position',[p_left+(.5*p_wid)+col*(p_wid+p_spacing) p_bot+row*(p_height+p_spacing) p_wid p_height])
            end
            var_plot_dirMax = nan_mean(dir_max_mean(:,:,i{1,x}),3)-(nan_mean(dir_max_mean(:,:,years1),3));
%             disp([min(var_plot_dir(:)) max(var_plot_dir(:))])
            pcolor(lon,lat,var_plot_dirMax) % plot
            shading interp
            set(gca,'FontSize', 8) % font size
            tit = sprintf('1950 - 2020 "Max" Dir vs %d "Max" Dir',starts(x));
            title(tit)
            caxis([-20 20])
            xlim([left_bound right_bound])
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
        chan.FontSize = 11;
        set(chan,'Position',[p_left .925 ((p_left + 4*p_wid + 3*p_spacing)-p_left) .02])
        printFig(gcf, 'MACA_Degree_Difference_Direction_CurMAX2Future', [11 8.5], 'png', 150)
    case 'avgdir'
        for x = 1:length(starts)
            if x <= 4
                axes('position',[p_left+col*(p_wid+p_spacing) p_bot+row*(p_height+p_spacing) p_wid p_height])
            else
                axes('position',[p_left+(.5*p_wid)+col*(p_wid+p_spacing) p_bot+row*(p_height+p_spacing) p_wid p_height])
            end
            var_plot_dir = nan_mean(dir_avg(:,:,i{1,x}),3)-(nan_mean(dir_avg(:,:,years1),3));
            disp([min(var_plot_dir(:)) max(var_plot_dir(:))]);
            pcolor(lon,lat,var_plot_dir) % plot
            shading flat
            set(gca,'FontSize', 8) % font size
            tit = sprintf('1950 - 2020 Avg Dir. vs %d Avg Dir.',starts(x));
            title(tit)
            caxis([-20 20])
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
        chan.FontSize = 11;
        set(chan,'Position',[p_left .925 ((p_left + 4*p_wid + 3*p_spacing)-p_left) .02])
        printFig(gcf, 'MACA_Degree_Difference_Direction_CurAVG2Future', [11 8.5], 'png', 150)
end