% This code does the following:
% 1. Loads in Obs data and finds MACA data closest to it
% 2. Make sure to save MACA and NNRP structures if you want to avoid
%   running this long finding process everytime
% 3. It plots the distributions of the data for overlapping times
% 4. Plots distributions of MACA and NNRP over overlapping historic time
% 5. Plots MACA station location average/max speed over 30 year segments
% 6. Plots distributions of MACA and Obs as well as future MACA 
% 7. Plots distributions of MACA at the Obs point in past and future 
addpath C:\Functions_Matlab
addpath C:\Functions_Matlab\time
clearvars
% Do you need to find the data closest to obs or have you already saved it?
data_find = 0;
switch data_find
    case 1
        % Load in station obs and grab lat/lon for station
        clearvars
        fol_loc = 'C:\Users\ahooshmand\Desktop\Data_Forcings\station_data\gap_hourly'; % location of data
        stations = {'bham_airport','mcchord_afb','vic_air','whidbey_nas'}; % station names
        for m = 1:length(stations)
            file_load = strcat(fol_loc,'\',stations{m},'_hourly.mat');
            load(file_load);
            obs(m).name = stations{m};
            obs(m).lat = lat;
            obs(m).lon = lon;
        end
        clear airtemp dewp file_load fol_loc lat lon m slp stp time wnddir wndspd stations
        
        
        % Load in MACA -----------------------------------
        tic
        % ---------------------- Historic Data ------------------------------------
        Mh_u = ncread('C:\Users\ahooshmand\Desktop\Data_Forcings\MACA\maca_u_GFDL_ESM2M_historical_1950_2005.nc','eastward_wind');
        Mh_v = ncread('C:\Users\ahooshmand\Desktop\Data_Forcings\MACA\maca_v_GFDL_ESM2M_historical_1950_2005.nc','northward_wind');
        Mh_time = ncread('C:\Users\ahooshmand\Desktop\Data_Forcings\MACA\maca_u_GFDL_ESM2M_historical_1950_2005.nc','time');
        Mh_time = double(Mh_time); % convert to double for following step
        Mh_time = datenum(1900,1,1+Mh_time); % Convert MACA time to Matlab Datenum
        Mh_lat = ncread('C:\Users\ahooshmand\Desktop\Data_Forcings\MACA\maca_u_GFDL_ESM2M_historical_1950_2005.nc','lat');
        Mh_lon = ncread('C:\Users\ahooshmand\Desktop\Data_Forcings\MACA\maca_u_GFDL_ESM2M_historical_1950_2005.nc','lon');
        Mh_lon = -1*(360 - Mh_lon);
        Mh_spd = hypot(Mh_u,Mh_v);
        wnddir_temp = (180/pi)*atan2(Mh_u,Mh_v);
        % Rotate winds to compass directions
        wnddir_temp = 90 - wnddir_temp;
        wnddir_temp(wnddir_temp<0)=wnddir_temp(wnddir_temp<0)+360;
        % switch to conventional coming from dir
        wnddir_temp = wnddir_temp+180;
        Mh_dir = wrapTo360(wnddir_temp);
        clear wnddir_temp
        
        % ---------------------- Future Data --------------------------------------
        scenario = 'rcp45';
        switch scenario
            case 'rcp45'
                Mf_u = ncread('C:\Users\ahooshmand\Desktop\Data_Forcings\MACA\maca_u_GFDL_ESM2M_rcp45_2006_2100.nc','eastward_wind');
                Mf_v = ncread('C:\Users\ahooshmand\Desktop\Data_Forcings\MACA\maca_v_GFDL_ESM2M_rcp45_2006_2100.nc','northward_wind');
                Mf_time = ncread('C:\Users\ahooshmand\Desktop\Data_Forcings\MACA\maca_u_GFDL_ESM2M_rcp45_2006_2100.nc','time');
                Mf_time = double(Mf_time); % convert to double for following step
                Mf_time = datenum(1900,1,1+Mf_time); % Convert MACA time to Matlab Datenum
                Mf_spd = hypot(Mf_u,Mf_v);
                wnddir_temp = (180/pi)*atan2(Mf_u,Mf_v);
                % Rotate winds to compass directions
                wnddir_temp = 90 - wnddir_temp;
                wnddir_temp(wnddir_temp<0)=wnddir_temp(wnddir_temp<0)+360;
                % switch to conventional coming from dir
                wnddir_temp = wnddir_temp+180;
                Mf_dir = wrapTo360(wnddir_temp);
                clear wnddir_temp
            case 'rcp85'
                Mf_u = ncread('C:\Users\ahooshmand\Desktop\Data_Forcings\MACA\maca_u_GFDL_ESM2M_rcp45_2006_2100.nc','eastward_wind');
                Mf_v = ncread('C:\Users\ahooshmand\Desktop\Data_Forcings\MACA\maca_v_GFDL_ESM2M_rcp45_2006_2100.nc','northward_wind');
                Mf_time = ncread('C:\Users\ahooshmand\Desktop\Data_Forcings\MACA\maca_u_GFDL_ESM2M_rcp45_2006_2100.nc','time');
                Mf_time = double(Mf_time); % convert to double for following step
                Mf_time = datenum(1900,1,1+Mf_time); % Convert MACA time to Matlab Datenum
                Mf_spd = hypot(Mf_u,Mf_v);
                wnddir_temp = (180/pi)*atan2(Mf_u,Mf_v);
                % Rotate winds to compass directions
                wnddir_temp = 90 - wnddir_temp;
                wnddir_temp(wnddir_temp<0)=wnddir_temp(wnddir_temp<0)+360;
                % switch to conventional coming from dir
                wnddir_temp = wnddir_temp+180;
                Mf_dir = wrapTo360(wnddir_temp);
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
        clear Mh_dir Mf_dir scenario
        toc
        
        % Find closest MACA point to obs station
        for x = 1:length(obs)
            dist_list = (M.lon-obs(x).lon).^2 + (M.lat-obs(x).lat).^2;
            [~,I] = min(dist_list(:));
            temp_x = M.lon(I);
            temp_y = M.lat(I);
            [r,c] = find(M.lon == temp_x & M.lat == temp_y);
            MA(x).u = M.u(r,c,:);
            MA(x).v = M.v(r,c,:);
            MA(x).time = M.time;
            MA(x).lat = M.lat(r,c);
            MA(x).lon = M.lon(r,c);
            MA(x).spd = M.spd(r,c,:);
            MA(x).dir = M.dir(r,c,:);
        end
        return
        
        % load NNRP -----------------------------------------
        tic
        N = load('C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\Storm_Tracker\west_coast_NNRP.mat');
        toc
        
        % Find closest NNRP point to obs station
        for x = 1:length(obs)
            dist_list = (N.lon_m-obs(x).lon).^2 + (N.lat_m-obs(x).lat).^2;
            [~,I] = min(dist_list(:));
            temp_x = N.lon_m(I);
            temp_y = N.lat_m(I);
            [r,c] = find(N.lon_m == temp_x & N.lat_m == temp_y);
            NN(x).time = N.time;
            NN(x).lat = N.lat_m(r,c);
            NN(x).lon = N.lon_m(r,c);
            NN(x).spd = N.wndspd(r,c,:);
            NN(x).dir = N.wnddir(r,c,:);
            NN(x).slp = N.slp(r,c,:);
        end
    case 0
        load('MACA_Obs_Points.mat')
        load('NNRP_Obs_Points.mat')
end
% Load in Obs
fol_loc = 'C:\Users\ahooshmand\Desktop\Data_Forcings\station_data\nan_hourly'; % location of data
stations = {'bham_airport','mcchord_afb','vic_air','whidbey_nas'}; % station names
for m = 1:length(stations)
    file_load = strcat(fol_loc,'/',stations{m},'_hourly.mat');
    load(file_load);
    % Grab meteo data for window of time
    O(m).name = stations{m};
    O(m).wndspd = wndspd;
    O(m).wnddir = wnddir';
    O(m).slp = slp;
    O(m).time = time;
    O(m).lat = lat;
    O(m).lon = lon;
    clear airtemp slp time wnddir wndspd lat lon
end
% Load in WA
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
clear airtemp dewp file_load fol_loc lat lon m slp stp time wnddir wndspd stations j wa_coast
%% Plot distributions for the different weather data - sub daily 

%nbins = 15;
x_locs = 0:0.5:25;
% xbins = vecotrize
% Reshape obs to be multi-row x 1-column matrix of daily averages - Note
% column number is the station number from obs list
for s = 1:length(O)
    temp_spd = O(s).wndspd;
    obs_day_vec = datenum(year(O(s).time(1)),month(O(s).time(1)),...
        day(O(s).time(1))):1:datenum(year(O(s).time(end)),...
        month(O(s).time(end)),day(O(s).time(end)));
    obs_temp_avg = zeros(length(obs_day_vec),1); % Initialize
    
    count = 0;
    for x = 1:round(length(temp_spd)/24)
        if x == 1
            i_st = 1;
            i_end = 24;
        else
            i_st = i_end +1;
            i_end = 24*x;
        end
        i_use = i_st:i_end;
        obs_temp_avg(x) = nan_mean(temp_spd(i_use));
    end
    
    % Find where MACA is same as obs
    ma_inds = findnearest(MA(s).time,obs_day_vec(1)):findnearest(MA(s).time,obs_day_vec(end));
    nn_inds = findnearest(NN(s).time,obs_day_vec(1)):findnearest(NN(s).time,obs_day_vec(end));
    
    % plot Distributions and save
    clf
    h1 = histogram(obs_temp_avg,x_locs,'Normalization','probability');
    l1 = h1.Values;
    hold on
    temp_ma = MA(s).spd(:,:,ma_inds);
    temp_ma = temp_ma(:);
    h2 = histogram(temp_ma,x_locs,'Normalization','probability','FaceAlpha',0.4);
    l2 = h2.Values;
    hold on
    h3 = histogram(NN(s).spd(:,:,nn_inds),x_locs,'Normalization','probability','FaceAlpha',0.4);
    l3 = h3.Values;
    legend('Obs','Maca','NNRP')
    xlabel('Wind Speed [m/s]')
    ylabel('Probability')
    outname = sprintf('MACA_vs_Obs_vs_NNRP_Hist_%s',O(s).name);
    printFig(gcf, outname, [14 14], 'png', 150)
    clf
    x_vec = 1:1:length(x_locs)-1;
    plot(x_vec,l1)
    hold on
    plot(x_vec,l2)
    hold on
    plot(x_vec,l3)
    legend('Obs','Maca','NNRP')
    xlabel('Wind Speed [m/s]')
    ylabel('Probability')
    fprintf('%s',O(s).name)
    outname = sprintf('MACA_vs_Obs_vs_NNRP_HistLINE_%s',O(s).name);
    printFig(gcf, outname, [14 14], 'png', 150)
    pause
end

%% Plot MACA NNRP and OBS on a Map to show proximity to one another
clf
p1 = plot(wa_lon, wa_lat, 'Color','k');
hold on
for x = 1:length(O)
    p2 = plot(O(x).lon,O(x).lat,'o','Color','k','MarkerFaceColor','k');
    hold on
    p3 = plot(NN(x).lon,NN(x).lat,'o','Color','r','MarkerFaceColor','r');
    hold on
    p4 = plot(MA(x).lon,MA(x).lat,'o','Color','b','MarkerFaceColor','b');
end
axis equal
legend([p2 p3 p4],'Obs','NNRP','MACA')
ylabel('Degrees Latitude')
xlabel('Degrees Longitude')
% Manually Zoom to see how close points are
%% Put NNRP and Obs on same timevec as MACA and then compare distributions
x_locs = 0:0.5:25;
for s = 1:length(O)
    % Convert Obs to Daily average
    temp_spd = O(s).wndspd;
    obs_day_vec = datenum(year(O(s).time(1)),month(O(s).time(1)),...
        day(O(s).time(1))):1:datenum(year(O(s).time(end)),...
        month(O(s).time(end)),day(O(s).time(end)));
    obs_temp_avg = zeros(length(obs_day_vec),1); % Initialize
    for x = 1:round(length(temp_spd)/24) % Find Indices and average
        if x == 1
            i_st = 1;
            i_end = 24;
        else
            i_st = i_end+1;
            i_end = 24*x;
        end
        i_use = i_st:i_end;
        obs_temp_avg(x) = nan_mean(temp_spd(i_use));
    end
    
    % Convert NNRP to Daily Average
    temp_spd = NN(s).spd(:,:,:);
    nn_day_vec = datenum(year(NN(s).time(1)),month(NN(s).time(1)),...
        day(NN(s).time(1))):1:datenum(year(NN(s).time(end)),...
        month(NN(s).time(end)),day(NN(s).time(end)));
    nn_temp_avg = zeros(length(nn_day_vec),1); % Initialize
    for x = 1:round(length(temp_spd)/4) % Find Indices and average
        if x == 1
            i_st = 1;
            i_end = 4;
        else
            i_st = i_end+1;
            i_end = 4*x;
        end
        i_use = i_st:i_end;
        nn_temp_avg(x) = nan_mean(temp_spd(:,:,i_use));
    end
    
    % Find where MACA is same as obs
    ma_inds = findnearest(MA(s).time,obs_day_vec(1)):findnearest(MA(s).time,obs_day_vec(end));
    nn_inds = findnearest(nn_day_vec,obs_day_vec(1)):findnearest(nn_day_vec,obs_day_vec(end));
    
    % plot Distributions and save
    clf
    h1 = histogram(obs_temp_avg,x_locs,'Normalization','probability');
    l1 = h1.Values;
    hold on
    temp_ma = MA(s).spd(:,:,ma_inds);
    temp_ma = temp_ma(:);
    h2 = histogram(temp_ma,x_locs,'Normalization','probability','FaceAlpha',0.4);
    l2 = h2.Values;
    hold on
    h3 = histogram(nn_temp_avg,x_locs,'Normalization','probability','FaceAlpha',0.4);
    l3 = h3.Values;
    legend('Obs - Daily Avg','Maca','NNRP - Daily Avg')
    xlabel('Wind Speed [m/s]')
    ylabel('Probability')
    outname = sprintf('MACA_vs_Obs_vs_NNRP_Daily_AVG_Hist_%s',O(s).name);
    printFig(gcf, outname, [14 14], 'png', 150)
    clf
    x_vec = 1:1:length(x_locs)-1;
    plot(x_vec,l1)
    hold on
    plot(x_vec,l2)
    hold on
    plot(x_vec,l3)
    legend('Obs - Daily Avg','Maca','NNRP - Daily Avg')
    xlabel('Wind Speed [m/s]')
    ylabel('Probability')
    fprintf('%s\n',O(s).name)
    outname = sprintf('MACA_vs_Obs_vs_NNRP_Daily_AVG_HistLINE_%s',O(s).name);
    printFig(gcf, outname, [14 14], 'png', 150)
    pause
end

%% Plot trends in stations vs MACA vs NNRP

% Find Climate Filters (30 year filters)
% MACA
starts1 = 1950:10:2070;
ends1 = 1980:10:2100;
im = cell(length(MA),length(starts1));
for xx = 1:length(starts1)
    for x = 1:length(MA)
        im{x,xx} = find(year(MA(x).time) >= starts1(xx) & year(MA(x).time) <= ends1(xx));
    end
end

% NNRP & NNRP
starts2 = 1950:10:1980;
ends2 = 1980:10:2010;
in = cell(length(MA),length(starts2));
for xx = 1:length(starts2)
    for x = 1:length(NN)
        in{x,xx} = find(year(NN(x).time) >= starts2(xx) & year(NN(x).time) <= ends2(xx));
        io{x,xx} = find(year(O(x).time) >= starts2(xx) & year(O(x).time) <= ends2(xx));
    end
end

% Find averages over Data
stns = [1,2,4]; % skip out on Victoria Airport
for x = 1:length(stns)
    % MACA First
    for j = 1:length(starts1)
        temp_max = max(MA(stns(x)).spd(:,:,im{x,j}),[],3);
        temp_mean = nan_mean(MA(stns(x)).spd(:,:,im{x,j}),3);
        m_max{j} = max(temp_max(:));
        m_mean{j} = nan_mean(temp_mean(:));
    end
    for j = 1:length(starts2)
        temp_nmax = max(NN(stns(x)).spd(:,:,in{x,j}),[],3);
        n_max{j} = max(temp_nmax(:));
        temp_omax = max(O(stns(x)).wndspd(io{x,j}));
        o_max{j} = max(temp_omax(:));
        temp_nmean = nan_mean(NN(stns(x)).spd(:,:,in{x,j}),3);
        n_mean{j} = nan_mean(temp_nmean(:));
        temp_omean = nan_mean(O(stns(x)).wndspd(io{x,j}));
        o_mean{j} = nan_mean(temp_omean(:));
    end
    
    
    % Plot
    clf
    t = 1:13;
    t2 = 1:4;
    metric = 'mean';
    switch metric
        case 'max'
            for x = 1:13
                ma_max(x) = m_max{1,x};
            end
            plot(t,ma_max)
            hold on
            for x = 1:4
                nn_max(x) = n_max{1,x};
                ob_max(x) = o_max{1,x};
            end
            plot(t2,nn_max)
            hold on
            plot(t2,ob_max)
            
            set(gca,'xtick',1:13)
            xl = {'1950-1980','1960-1990','1970-2000','1980-2010','1990-2020',...
                '2000-2030','2010-2040','2020-2050','2030-2060','2040-2070','2050-2080','2060-2090','2070-2100'};
            set(gca,'xticklabel',xl)
            set(gca,'fontsize',8)
            xlabel('Time - 30 year segments')
            ylabel('Max Windspeed [m/s]')
            grid on
            legend('MACA','NNRP','OBS')
            pause
        case 'mean'
            for x = 1:13
                ma_mean(x) = m_mean{1,x};
            end
            plot(t,ma_mean)
            hold on
            for x = 1:4
                nn_mean(x) = n_mean{1,x};
                ob_mean(x) = o_mean{1,x};
            end
            plot(t2,nn_mean)
            hold on
            plot(t2,ob_mean)
            
            set(gca,'xtick',1:13)
            xl = {'1950-1980','1960-1990','1970-2000','1980-2010','1990-2020',...
                '2000-2030','2010-2040','2020-2050','2030-2060','2040-2070','2050-2080','2060-2090','2070-2100'};
            set(gca,'xticklabel',xl)
            set(gca,'fontsize',8)
            xlabel('Time - 30 year segments')
            ylabel('Max Windspeed [m/s]')
            grid on
            legend('MACA','NNRP','OBS')
            pause
    end
end
%% 

x_locs = 0:0.5:25;
for s = 1:length(O)
    % Convert Obs to Daily average
    temp_spd = O(s).wndspd;
    obs_day_vec = datenum(year(O(s).time(1)),month(O(s).time(1)),...
        day(O(s).time(1))):1:datenum(year(O(s).time(end)),...
        month(O(s).time(end)),day(O(s).time(end)));
    obs_temp_avg = zeros(length(obs_day_vec),1); % Initialize
    for x = 1:round(length(temp_spd)/24) % Find Indices and average
        if x == 1
            i_st = 1;
            i_end = 24;
        else
            i_st = i_end+1;
            i_end = 24*x;
        end
        i_use = i_st:i_end;
        obs_temp_avg(x) = nan_mean(temp_spd(i_use));
    end
    
    % Convert NNRP to Daily Average
    temp_spd = NN(s).spd(:,:,:);
    nn_day_vec = datenum(year(NN(s).time(1)),month(NN(s).time(1)),...
        day(NN(s).time(1))):1:datenum(year(NN(s).time(end)),...
        month(NN(s).time(end)),day(NN(s).time(end)));
    nn_temp_avg = zeros(length(nn_day_vec),1); % Initialize
    for x = 1:round(length(temp_spd)/4) % Find Indices and average
        if x == 1
            i_st = 1;
            i_end = 4;
        else
            i_st = i_end+1;
            i_end = 4*x;
        end
        i_use = i_st:i_end;
        nn_temp_avg(x) = nan_mean(temp_spd(:,:,i_use));
    end
    
    % Find where MACA is same as obs
    ma_inds = findnearest(MA(s).time,obs_day_vec(1)):findnearest(MA(s).time,obs_day_vec(end));
    nn_inds = findnearest(nn_day_vec,obs_day_vec(1)):findnearest(nn_day_vec,obs_day_vec(end));
    fut_inds = find(year(MA(s).time) >= 2020 & year(MA(s).time) <= 2100);
    % plot Distributions and save
    clf
    h1 = histogram(obs_temp_avg,x_locs,'Normalization','probability');
    l1 = h1.Values;
    hold on
    temp_ma = MA(s).spd(:,:,ma_inds);
    temp_ma = temp_ma(:);
    h2 = histogram(temp_ma,x_locs,'Normalization','probability','FaceAlpha',0.4);
    l2 = h2.Values;
    hold on
    h3 = histogram(nn_temp_avg,x_locs,'Normalization','probability','FaceAlpha',0.4);
    l3 = h3.Values;
    temp_ma = MA(s).spd(:,:,fut_inds);
    temp_ma = temp_ma(:);
    h4 = histogram(temp_ma,x_locs,'Normalization','probability','FaceAlpha',0.4);
    l4 = h4.Values;
    legend('Obs - Daily Avg','Maca','NNRP - Daily Avg','Future Maca')
    xlabel('Wind Speed [m/s]')
    ylabel('Probability')
    %outname = sprintf('MACA_vs_Obs_vs_NNRP_Daily_AVG_Hist_%s',O(s).name);
    %printFig(gcf, outname, [14 14], 'png', 150)
    clf
    x_vec = 1:1:length(x_locs)-1;
    plot(x_vec,l1)
    hold on
    plot(x_vec,l2)
    hold on
    plot(x_vec,l3)
    hold on 
    plot(x_vec,l4)
    legend('Obs - Daily Avg','Maca','NNRP - Daily Avg','Future MACA')
    xlabel('Wind Speed [m/s]')
    ylabel('Probability')
    fprintf('%s\n',O(s).name)
    outname = sprintf('%s_Hist_Compare_FutureMACA',O(s).name);
    printFig(gcf, outname, [14 14], 'png', 150)
    pause
end

%% Plot distribution of MACA at Obs spots in past and future

stns = [1 2 4];
xlocs = 0:1:15;
for x = 1:length(stns)
    past_inds = find(year(MA(stns(x)).time) >= 1950 & year(MA(x).time) <= 2020);
    fut_inds = find(year(MA(stns(x)).time) >= 2020 & year(MA(x).time) <= 2100);
%     spd_past = MA.spd(:,:,past_inds); spd_past = spd_past(:);
%     spd_fut = MA.spd(:,:,fut_inds); spd_fut = spd_fut(:);
    clf
    h1 = histogram(MA(stns(x)).spd(:,:,past_inds),x_locs,'Normalization','probability');
    hold on
    h2 = histogram(MA(stns(x)).spd(:,:,fut_inds),x_locs,'Normalization','probability','FaceAlpha',0.4);
    legend('MACA Historic','MACA Future')
    xlabel('Wind Speed [m/s]')
    ylabel('Probability')
    outname = sprintf('MACA_at_%s_Past_vs_Future',O(stns(x)).name);
    printFig(gcf,outname,[8.5 11],'png',150)    
end

%% Print out the Percentiles 

for x = 1:length(stns)
    fprintf('Past 90th: %4.2f\n',prctile(MA(stns(x)).spd(:,:,past_inds),90));
    fprintf('Future 90th: %4.2f\n',prctile(MA(stns(x)).spd(:,:,fut_inds),90));
    fprintf('Past 99th: %4.2f\n',prctile(MA(stns(x)).spd(:,:,past_inds),99));
    fprintf('Future 99th: %4.2f\n',prctile(MA(stns(x)).spd(:,:,fut_inds),99));
    fprintf('\n')
    pause
end