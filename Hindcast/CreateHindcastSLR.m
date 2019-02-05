%   Creates a hindcast for a specific transect
%   This script will create a hindcast along any specific transect created
%   in google earth.  The user determines what type of data they want to
%   use as a forcing.  The data can be obs, quantile corrected NNRP, or
%   NNRP.  This function should be called in a script and then the output
%   should be used to plot.

clearvars
addpath C:\Functions_Matlab
slr_val = 3.0;
ft2m = .3048;
slr_val = slr_val*ft2m; %[m]
slr = 1;
if slr
    
    if rem((slr_val/ft2m),round(slr_val/ft2m)) == 0
        slrFolder = sprintf('%d_FT_SLR',slr_val/ft2m);
    else
        slrFolder = sprintf('%d_FT_SLR',100*(slr_val/ft2m));
    end
    if ~exist(sprintf('Hindcast_SLR/%s',slrFolder),'dir')
        mkdir(sprintf('Hindcast_SLR/%s',slrFolder),'dir')
    end
end

% Type of output that we want - Bulk vs TS
output = 'Bulk'; % Timeseries


% ----------------- Folder Locations of data ------------------------------
nnrp_fol = 'E:\Abbas\PS_COSMOS\Thesis_Modeling\Quantile_Correction\NNRP_Updated\';
mask_fol = 'E:\Abbas\PS_COSMOS\Thesis_Modeling\KML\Basin_Masks\';
tran_fol = 'E:\Abbas\PS_COSMOS\Thesis_Modeling\KML\PugetSoundShorline\';
lut_fol = 'E:\Abbas\PS_COSMOS\Thesis_Modeling\LUT\Offshore_LUT_Extracted\';
tide_fol = 'E:\Abbas\Model_Met_Forcings\Crosby_TWLReconstructions\Reconstructions\';

% ----------------- Key Value Pairs for Basin and NNRP --------------------
keys = {'Anacortes.kml','BellinghamBay.kml','Blaine.kml','Hood_Hook.kml',...
    'Hood_Kitsap.kml','Middle_Hood.kml','Port_Ludlow.kml','Nor_Hood.kml',...
    'LummiBay.kml','LummiIsland.kml','NorthCentral.kml','Olympia.kml','OlympiaNorth.kml','OlympiaWest.kml',...
    'PadillaBay.kml','PortTownsend_Hood.kml','Skagit.kml','SouthCentral.kml',...
    'Tacoma.kml','Whidbey_JDF.kml','WhidbeyAnacortes.kml','WhidbeyBasin.kml'...
    'JDF1.kml','JDF2.kml','JDF3.kml','JDF4.kml','JDF5.kml','JDF6.kml','JDF7.kml',...
    'JDF8.kml','JDF9.kml','JDF10.kml','JDF11.kml','GuemesCypress_Mask1.kml',...
    'CypressSinclairGuemes_Mask1.kml','Guemes_Mask1.kml','LopezIsland_Mask1.kml',...
    'Orcas_Mask1.kml','SanJuan_Mask1.kml','SanJuan_Mask2.kml','SanJuan_mask3.kml',...
    'Stuart_Mask1.kml','Waldron_Mask1.kml','Waldron_Mask2.kml','Sucia_Mask1.kml',...
    'BoundaryBay_Mask.kml','PointRoberts_Mask.kml','Tsawwassen_Mask.kml',...
    'RichmondBC_Mask.kml','IonaBeachBC_Mask.kml','BurrardInlet_Mask.kml',...
    'BowenIslandBC_Mask.kml','SmithIsland.kml'};

values = {'Cypress_Guemes_NNRP_Point.mat','Chuckanut_NNRP_Point.mat',...
    'Blaine_NNRP_Point.mat','Hood_Hook_CORRECTED_NNRP_Point.mat',...
    'Hood_Kitsap_CORRECTED_NNRP_Point.mat','Middle_Hood_Canal_CORRECTED_NNRP_Point.mat',...
    'Port_Ludlow_CORRECTED_NNRP_Point.mat','Nor_Hood_Canal_CORRECTED_NNRP_Point.mat',...
    'Lummi_Bay_NNRP_Point.mat','Lummi_Island_NNRP_Point.mat','Kingston_NNRP_Point.mat',...
    'McNeil_Island_NNRP_Point.mat','Nor_Olympia_CORRECTED_NNRP_Point.mat',...
    'Middle_Olympia_CORRECTED_NNRP_Point.mat','Padilla_Bay_NNRP_Point.mat',...
    'Townsend_Whidbey_NNRP_Point.mat','Skagit_Bay_NNRP_Point.mat','Blake_Island_NNRP_Point.mat',...
    'Tacoma_CORRECTED_NNRP_Point.mat','JDF_Whidbey_NNRP_Point.mat','Deception_Pass_NNRP_Point.mat',...
    'Whidbey_Camano_NNRP_Point.mat','JDF1_NNRP_Point.mat','Neah_Bay_NNRP_Point.mat',...
    'JDF3_NNRP_Point.mat','JDF4_NNRP_Point.mat','JDF5_NNRP_Point.mat','JDF6_NNRP_Point.mat',...
    'JDF7_NNRP_Point.mat','JDF8_NNRP_Point.mat','JDF9_NNRP_Point.mat','Dungeness_NNRP_Point.mat',...
    'JDF11_NNRP_Point.mat','Cypress_Guemes_NNRP_Point.mat','Lummi_Island_NNRP_Point.mat',...
    'Padilla_Bay_NNRP_Point.mat','Lopez_South_NNRP_Point.mat','SJ_Islands_NNRP_Point.mat',...
    'SanJuan_South_NNRP_Point.mat','SanJuan_Island_NNRP_Point.mat','Stuart_Island_NNRP_Point.mat',...
    'SaturnaIsland_NNRP_Point.mat','Orcas_Island_NNRP_Point.mat','Saturna_Waldron_NNRP_Point.mat',...
    'Sucia_Island_NNRP_Point.mat','Boundary_Bay_NNRP_Point.mat','Point_Roberts_NNRP_Point.mat',...
    'Tsawwassen_NNRP_Point.mat','RichmondBC_NNRP_Point.mat','IonaBeach_BC_NNRP_Point.mat',...
    'BurrardInlet_BC_NNRP_Point.mat','BowenIsland_BC_NNRP_Point.mat','Smith_Island_NNRP_Point.mat'};
% Sort the keys alphabetically
[keys,I] = sort(keys);
values = values(I);

DM = containers.Map(keys,values); % Creates a dictionary of key value pairs

% ------------------ Load LUT ---------------------------------------------
% L1 = load([lut_fol 'SWAN_10m_JDFLUT_offshore_extract.mat']);
% L2 = load([lut_fol 'SWAN_10m_PSLUT_offshore_extract.mat']);
% L3 = load([lut_fol 'SWAN_10m_SGALUT_offshore_extract.mat']);
LUT_list = dir([lut_fol '*.mat']);
for ii = 1:length(LUT_list)
    J(ii) = load([lut_fol LUT_list(ii).name]);
end
for ii = 1:length(LUT_list)
    J(ii).name = LUT_list(ii).name;
end
% ----------------- Load in single NNRP Point for the time vector ---------
nn2load = [nnrp_fol 'McNeil_Island_NNRP_Point.mat'];
N = load(nn2load);
N.time  = N.time';
n_time = N.time;
% Interp NNRP from 6 hourly to hourly
time = datenum(n_time(1):1/24:n_time(end));
clear N

% % Find each year and house those indices in a variable
years = year(time(1)):1:year(time(end));
yrInds = cell(length(years),1);
for ii = 1:length(yrInds)
    yrInds{ii} = find(year(time) == years(ii));
end
%
% % Find each month in the timeseries and house that in a variable
% moInds = cell(12,1);
% for ii = 1:length(moInds)
%     moInds{ii} = find(month(time) == ii);
% end


% ------------------ Load in Tides ----------------------------------------
tide_stns = dir([tide_fol '*.mat']);
for ii = 1:length(tide_stns)
    if strcmp(tide_stns(ii).name,'seattle_hrV_NAVD88.mat')
        load([tide_fol tide_stns(ii).name]);
        T(ii).time = tides.time;
        T(ii).twl = double(tides.WL_VALUE);
        T(ii).lon = tides.longitude(1);
        T(ii).lat = tides.latitude(1);
        T(ii).name = tide_stns(ii).name;
        
    else
        temp = load([tide_fol tide_stns(ii).name]);
        T(ii).time = temp.time;
        T(ii).twl = temp.twl;
        T(ii).name = tide_stns(ii).name;
        T(ii).lon = temp.lon;
        T(ii).lat = temp.lat;
    end
    % Interpolate twl onto time vec
    inds = find(unique(T(ii).time));
    T(ii).time = T(ii).time(inds);
    T(ii).twl = T(ii).twl(inds);
    inds = find(diff(T(ii).time) == 0);
    T(ii).time(inds) = [];
    T(ii).twl(inds) = [];
    T(ii).twl_i = double(interp1(T(ii).time,T(ii).twl,time));
end

% ----------------- Load Transect -----------------------------------------
% K = kml2struct([tran_fol 'Swan_10mContour.kml']);

% ----------------- Load in Masks -----------------------------------------
list = dir('E:\Abbas\PS_COSMOS\Thesis_Modeling\KML\Basin_Masks\*.kml');
in = [];
for ii = 1:length(list)
    M(ii) = kml2struct([mask_fol list(ii).name]); % Load in the basin KML
    M(ii).Lon(end) = []; M(ii).Lat(end) = [];
end
% Convert Masks to UTM - not sure why it wouldn't let me put it in the loop
% but oh well it works

for ii = 1:length(M)
    [M(ii).x,M(ii).y] = deg2utm(M(ii).Lat,M(ii).Lon);
end
% Hout = dir('Hindcast_Output\*.mat');
disp('Done Loading Data - Beginning Hindcast')

% Load tide masks
maskFol = 'E:\Abbas\PS_COSMOS\Thesis_Modeling\KML\TideStationMasks\';
f = dir([maskFol '*.kml']);
for ii = 1:length(f)
    temp = kml2struct([maskFol f(ii).name]);
    K(ii).name = f(ii).name;
    K(ii).lon = temp.Lon;
    K(ii).lat = temp.Lat;
    [K(ii).x,K(ii).y] = deg2utm(temp.Lat,temp.Lon);
end

%% Create Hindcast from LUT
% ------------- Dimensions of LUT -----------------------------------------
% ------------- setup for interp ------------------------------------------
[X,Y,Z] = meshgrid(J(1).tide,J(1).speed,J(1).direc);
n_inds = [];
% ------------- Loop over basins then extraction points -------------------
tic
for kk = 1:length(keys) % Loop through the basin masks and grab subset of points to make hindcast of
    clear L
    [~,basinName,~] = fileparts(keys{kk});
    
    % See What parts of each LUT are within the Basin Mask
    in = cell(length(J),1);
    for ii = 1:length(J)
        in{ii} = inpolygon(J(ii).lon,J(ii).lat,M(kk).Lon,M(kk).Lat);
    end
    inList = cellfun(@sum,in); % get total points in each LUT in the current mask
    inInds = find(inList > 0); % grab indices of LUTs in mask
    if length(inInds) == 1 % Grab the LUT to use if it's only 1
        %         L = J(inInds);
        temp = J(inInds);
        tranInds = in{inInds}; % Subsample transect to be only part in basin
        L.hs = temp.hs(tranInds,:,:,:);
        L.tm = temp.tm(tranInds,:,:,:);
        L.tp = temp.tp(tranInds,:,:,:);
        L.x = temp.x(tranInds);
        L.y = temp.y(tranInds);
        L.lat = temp.lat(tranInds);
        L.lon = temp.lon(tranInds);
        L.hdir = temp.hdir(tranInds,:,:,:);
    else
        % Otherwise Find which parts of each LUT are in which basin and combine into a
        % single structure
        for jj = 1:length(inInds)
            temp = J(inInds(jj));
            tranInds = in{inInds(jj)}; % Subsample transect to be only part in basin
            t(jj).hs = temp.hs(tranInds,:,:,:);
            t(jj).tm = temp.tm(tranInds,:,:,:);
            t(jj).tp = temp.tp(tranInds,:,:,:);
            t(jj).x = temp.x(tranInds);
            t(jj).y = temp.y(tranInds);
            t(jj).lat = temp.lat(tranInds);
            t(jj).lon = temp.lon(tranInds);
            t(jj).hdir = temp.hdir(tranInds,:,:,:);
        end
        % Combine in to one structure
        L = struct;
        L.hs = [];
        L.tm = [];
        L.tp = [];
        L.x = [];
        L.y = [];
        L.lon = [];
        L.lat = [];
        L.hdir = [];
        for ll = 1:length(t)
            L.hs = cat(1,L.hs,t(ll).hs);
            %             L.hs = cat(L.hs,t(ll).hs);
            L.tm = cat(1,L.tm,t(ll).tm);
            L.tp = cat(1,L.tp,t(ll).tp);
            L.x = vertcat(L.x,t(ll).x);
            L.y = vertcat(L.y,t(ll).y);
            L.lon = vertcat(L.lon,t(ll).lon);
            L.lat = vertcat(L.lat,t(ll).lat);
            L.hdir = vertcat(L.hdir,t(ll).hdir);
        end
    end
    
    % ------------------- Load Forcings -----------------------------------
    % Grab NNRP Forcing
    NN_point = DM(keys{kk});
    N = load([nnrp_fol NN_point]);
    %     N.wnddir = squeeze(N.wnddir);
    %     N.wndspd = squeeze(N.wndspd);
    
    % Convert to U and V space -  Already done for me
    %     [N.U,N.V] = wind2UV('from',N.wnddir,'compass',N.wndspd);
    U = interp1(N.time,N.wind_u,time)';
    V = interp1(N.time,N.wind_v,time)';
    % Convert Back to wnddir and speed
    %     speed = sqrt(U.^2+V.^2);
    %     wnddir = wrap2360(270-(atan2(V,U)*(180/pi)));
    [speed,wnddir] = uv2compass(U,V);
    
    
    % --------- Loop through the LUT and interp to make Hindcast ----------
    for ff = 1:length(L.x)
        outName = sprintf('%s_%d_%d_waveOut.mat',basinName,ff,length(L.x));
        if ~exist(['E:\Abbas\PS_COSMOS\Thesis_Modeling\LUT\Hindcast_SLR\' slrFolder '\' outName])
            fprintf('%s doesnt exist, making Hindcast\n',basinName)
            
            % Find which tide station to use
            mask2use = zeros(length(K),1);
            for tt = 1:length(K)
                mask2use(tt) = inpolygon(L.x(ff),L.y(ff),K(tt).lon,K(tt).lat);
            end
            if sum(mask2use) == 1
                ping = find(mask2use ~= 0);
                % Set TWL based on ping
                switch ping
                    case 1 % Cherry
                        twl = T(3).twl_i';
                    case 2 % Friday
                        twl = T(4).twl_i';
                    case 3 % Neah
                        twl = T(5).twl_i';
                    case 4 % Port Angeles
                        twl = T(6).twl_i';
                    case 5 % Port Townsend
                        twl = T(1).twl_i';
                    case 6 % Seattle
                        twl = T(7).twl_i';
                    case 7 % Tacoma
                        twl = T(2).twl_i';
                end
                clear ping
                % If the point is not within a poly mask, find the closest tide
                % gauge and use that
            elseif sum(mask2use) == 0 || sum(mask2use) > 1
                dists = zeros(length(T),1);
                for ss = 1:length(T)
                    dists(ss) = sqrt((L.x(ff) - T(ss).lon).^2 + (L.y(ff) - T(ss).lat).^2);
                end
                [~,I] = min(dists);
                twl = T(I).twl_i';
            end
            
            
            % Initialize Variables
            hs_ts = zeros(1,length(time));
            tp_ts = hs_ts;
            dir_ts = hs_ts;
            tm_ts = hs_ts;
            
            temp = permute(squeeze(L.hs(ff,:,:,:)),[2 1 3]);
            nan_inds = isnan(temp);
            temp(nan_inds) = 0;
            if ~isnan(temp(2,1))
                if slr
                    hs_ts(1,:) = interp3(X,Y,Z,temp,twl+slr_val,speed,wnddir,'linear')';
                else
                    hs_ts(1,:) = interp3(X,Y,Z,temp,twl,speed,wnddir,'linear')'; % Hsig timeseries
                end
                %---Peak Period
                temp = permute(squeeze(L.tp(ff,:,:,:)),[2 1 3]);
                nan_inds = isnan(temp);
                temp(nan_inds) = 0;
                if slr
                    tp_ts(1,:) = interp3(X,Y,Z,temp,twl+slr_val,speed,wnddir,'linear')';
                else
                    tp_ts(1,:) = interp3(X,Y,Z,temp,twl,speed,wnddir,'linear')';
                end
                %---Mean Period
                temp = permute(squeeze(L.tm(ff,:,:,:)),[2 1 3]);
                nan_inds = isnan(temp);
                temp(nan_inds) = 0;
                if slr
                    tm_ts(1,:) = interp3(X,Y,Z,temp,twl+slr_val,speed,wnddir,'linear')';
                else
                    tm_ts(1,:) = interp3(X,Y,Z,temp,twl,speed,wnddir,'linear')';
                end
                %---Wave Direction
                temp = permute(squeeze(L.hdir(ff,:,:,:)),[2 1 3]);
                nan_inds = isnan(temp);
                temp(nan_inds) = 0;
                if slr
                    dir_ts(1,:) = interp3(X,Y,Z,temp,twl+slr_val,speed,wnddir,'linear')';
                else
                    dir_ts(1,:) = interp3(X,Y,Z,temp,twl,speed,wnddir,'linear')';
                end
                fprintf('Completed %d out of %d...Moving On\n',ff,length(L.x))
            else
                fprintf('%d out of %d is a NAN point...Moving On',ff,length(L.x))
            end
            % ------------ Set any negative hs to zero (caused by winds < 2.5m/s) -----
            hs_ts(hs_ts<0) = 0;
            
            % ------------ Set any negative hs to zero (caused by winds < 2.5m/s) -----
            tp_ts(tp_ts<0) = 0;
            
            % ------------ Set any negative hs to zero (caused by winds < 2.5m/s) -----
            dir_ts(dir_ts<0) = 0;
            
            tm_ts(tm_ts<0) = 0;
            
            % Create vectors of U and V from wave height and direction
            % Note that we are using energy in lue of wave height, E is
            % proportional to hs squared
            [u,v] = compass2uv(hs_ts.^2,dir_ts);
            
            switch output
                case 'Bulk'
                    % Stats to compute
                    % 1. Daily Max
                    % 2. Daily Mean
                    % 3. Avg. Annual Max
                    % 4. Avg. Annual Mean
                    % 5. Monthly Max
                    % 6. Monthly Mean
                    % 7. Max Hs
                    % 8. Mean Hs
                    
                    % Make a time table of the Hsig 
                    ttime = datetime(time,'ConvertFrom','datenum');
                    TT = timetable(ttime',hs_ts',tp_ts',tm_ts',u',v');
                    % Calculate Statistics about the data
%                     D.dailyMax = mean(table2array(retime(TT,'daily','max')));
%                     D.dailyAvg = mean(table2array(retime(TT,'daily','mean')));
                    D.annualMax = mean(table2array(retime(TT,'yearly','max')));
                    D.annualMean = mean(table2array(retime(TT,'yearly','mean')));
                    monthlyMax = retime(TT,'monthly','max');
                    monthlyMean = retime(TT,'monthly','mean');
                    % Find monthly indices
                    ja = 1:12:height(monthlyMax);
                    fb = ja+1;
                    mo = fb + 1;
                    ap = mo + 1;
                    ma = ap + 1;
                    ju = ma + 1;
                    jl = ju + 1;
                    au = jl + 1;
                    sp = au + 1;
                    oc = sp + 1;
                    no = oc + 1;
                    dc = no + 1;
                    
                    D.janMax = mean(table2array(monthlyMax(ja,:)));
                    D.janMean = mean(table2array(monthlyMean(ja,:)));
                    D.febMax = mean(table2array(monthlyMax(fb,:)));
                    D.febMean = mean(table2array(monthlyMean(fb,:)));
                    D.marMax = mean(table2array(monthlyMax(mo,:)));
                    D.marMean = mean(table2array(monthlyMean(mo,:)));
                    D.aprMax = mean(table2array(monthlyMax(ap,:)));
                    D.aprMean = mean(table2array(monthlyMean(ap,:)));
                    D.mayMax = mean(table2array(monthlyMax(ma,:)));
                    D.mayMean = mean(table2array(monthlyMean(ma,:)));
                    D.junMax = mean(table2array(monthlyMax(ju,:)));
                    D.junMean = mean(table2array(monthlyMean(ju,:)));
                    D.julMax = mean(table2array(monthlyMax(jl,:)));
                    D.julMean = mean(table2array(monthlyMean(jl,:)));
                    D.augMax = mean(table2array(monthlyMax(au,:)));
                    D.augMean = mean(table2array(monthlyMean(au,:)));
                    D.sepMax = mean(table2array(monthlyMax(sp,:)));
                    D.sepMean = mean(table2array(monthlyMean(sp,:)));
                    D.octMax = mean(table2array(monthlyMax(oc,:)));
                    D.octMean = mean(table2array(monthlyMean(oc,:)));
                    D.novMax = mean(table2array(monthlyMax(no,:)));
                    D.novMean = mean(table2array(monthlyMean(no,:)));
                    D.decMax = mean(table2array(monthlyMax(dc,:)));
                    D.decMean = mean(table2array(monthlyMean(dc,:)));

                    D.tsMax = nanmax(hs_ts);
                    D.tsMean = nanmean(hs_ts);
                    
                    % Check on stats above
                    %                 temp_yr = zeros(1,length(years));
                    %                 temp_yrmean = temp_yr;
                    %                 for yy = 1:length(yrInds)
                    %                     inds = yrInds{yy};
                    %                     hs_yr_temp = hs_ts(1,inds);
                    %                     [temp_yr(yy), ~] = max(hs_yr_temp,[],2);
                    %                     temp_yrmean(yy) = mean(hs_yr_temp);
                    %                 end
                    %                 temp_yr = mean(temp_yr,2);
                    %                 temp_mo = mean(temp_yrmean,2);
                    
                    % Save The Entire Timeseries
                    D.lon = L.lon(ff);
                    D.lat = L.lat(ff);
                    D.slrMeters = slr_val;
                    D.slrFeet = slr_val/ft2m;
                    D.DataFormat = 'Hsig, Peak Period, Mean Period, U, V';
                    
                    save(outName,'-Struct','D')
                    disp('Wave Hindcast Saved...Moving On')
                    
                    movefile(outName,['E:\Abbas\PS_COSMOS\Thesis_Modeling\LUT\Hindcast_SLR\' slrFolder])
                    
                case 'Timeseries'
                    % Save The Entire Timeseries
                    lon = L.lon(ff);
                    lat = L.lat(ff);
                    slrMeters = slr_val;
                    slrFeet = slr_val/ft2m;
                    
                    save(outName,'lon','lat','hs_ts','tp_ts','tm_ts','dir_ts','speed','wnddir',...
                        'twl','time','slrMeters','slrFeet')
                    disp('Wave Hindcast Saved...Moving On')
                    movefile(outName,'E:\Abbas\PS_COSMOS\Thesis_Modeling\LUT\Hindcast_Output')
            end
            
            % Status Message
            fprintf('Done with %s: Zone %d of %d\n',basinName,ff,length(L.x))
        end
    end
    fprintf('%s Fully Completed - %.1f Percent of Basins\n',basinName,round((kk/length(keys)*100),1))
end

toc
%% Old Code I don't want to delete just yet

% Find which LUT to use and if they're overlapping use specific part of both
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     if ~isempty(inInds1) && isempty(inInds2) && isempty(inInds3)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         L = L1;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         disp('Using LUT Zone 1')
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         tranInds = inInds1; % Subsample transect to be only part in basin
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         L.hs = L.hs(tranInds,:,:,:);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         L.tm = L.tm(tranInds,:,:,:);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         L.tp = L.tp(tranInds,:,:,:);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         L.x = L.x(tranInds);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         L.y = L.y(tranInds);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         L.x_u = L.x_u(tranInds);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         L.y_u = L.y_u(tranInds);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         L.lat = L.lat(tranInds);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         L.lon = L.lon(tranInds);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     elseif ~isempty(inInds1) && ~isempty(inInds2) && isempty(inInds3)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         % Subsample first transect to get only part in basin
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         t1 = L1;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         t1.hs = t1.hs(inInds1,:,:,:);t1.hs = flipud(t1.hs);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         t1.tm = t1.tm(inInds1,:,:,:);t1.tm = flipud(t1.tm);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         t1.tp = t1.tp(inInds1,:,:,:);t1.tp = flipud(t1.tp);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         t1.x = t1.x(inInds1);t1.x = flipud(t1.x);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         t1.y = t1.y(inInds1);t1.y = flipud(t1.y);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         t1.x_u = t1.x_u(inInds1);t1.x_u = flipud(t1.x_u);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         t1.y_u = t1.y_u(inInds1);t1.y_u = flipud(t1.y_u);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         t1.lat = t1.lat(inInds1);t1.lat = flipud(t1.lat);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         t1.lon = t1.lon(inInds1);t1.lon = flipud(t1.lon);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         % Subsample second transect
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         t2 = L2;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         t2.hs = t2.hs(inInds2,:,:,:);t2.hs = flipud(t2.hs);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         t2.tm = t2.tm(inInds2,:,:,:);t2.tm = flipud(t2.tm);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         t2.tp = t2.tp(inInds2,:,:,:);t2.tp = flipud(t2.tp);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         t2.x = t2.x(inInds2);t2.x = flipud(t2.x);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         t2.y = t2.y(inInds2);t2.y = flipud(t2.y);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         t2.x_u = t2.x_u(inInds2);t2.x_u = flipud(t2.x_u);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         t2.y_u = t2.y_u(inInds2);t2.y_u = flipud(t2.y_u);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         t2.lat = t2.lat(inInds2);t2.lat = flipud(t2.lat);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         t2.lon = t2.lon(inInds2);t2.lon = flipud(t2.lon);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         % Combine in to one structure
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         L = struct;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         L.hs = vertcat(t1.hs,t2.hs);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         L.tm = vertcat(t1.tm,t2.tm);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         L.tp = vertcat(t1.tp,t2.tp);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         L.x = vertcat(t1.x,t2.x);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         L.y = vertcat(t1.y,t2.y);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         L.x_u = vertcat(t1.x_u,t2.x_u);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         L.y_u = vertcat(t1.y_u,t2.y_u);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         L.lon = vertcat(t1.lon,t2.lon);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         L.lat = vertcat(t1.lat,t2.lat);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     elseif ~isempty(inInds2) && isempty(inInds1) && isempty(inInds3)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         L = L2;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         disp('Using LUT Zone 2')
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         tranInds = inInds2; % Subsample transect to be only part in basin
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         L.hs = L.hs(tranInds,:,:,:);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         L.tm = L.tm(tranInds,:,:,:);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         L.tp = L.tp(tranInds,:,:,:);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         L.x = L.x(tranInds);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         L.y = L.y(tranInds);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         L.x_u = L.x_u(tranInds);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         L.y_u = L.y_u(tranInds);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         L.lat = L.lat(tranInds);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         L.lon = L.lon(tranInds);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     elseif ~isempty(inInds2) && ~isempty(inInds3) && isempty(inInds1)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %       % Subsample first transect to get only part in basin
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %       t1 = L2;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %       t1.hs = t1.hs(inInds2,:,:,:);t1.hs = flipud(t1.hs);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %       t1.tm = t1.tm(inInds2,:,:,:);t1.tm = flipud(t1.tm);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %       t1.tp = t1.tp(inInds2,:,:,:);t1.tp = flipud(t1.tp);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %       t1.x = t1.x(inInds2);t1.x = flipud(t1.x);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %       t1.y = t1.y(inInds2);t1.y = flipud(t1.y);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %       t1.x_u = t1.x_u(inInds2);t1.x_u = flipud(t1.x_u);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %       t1.y_u = t1.y_u(inInds2);t1.y_u = flipud(t1.y_u);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %       t1.lat = t1.lat(inInds2);t1.lat = flipud(t1.lat);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %       t1.lon = t1.lon(inInds2);t1.lon = flipud(t1.lon);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %       % Subsample second transect
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %       t2 = L3;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %       t2.hs = t2.hs(inInds3,:,:,:);t2.hs = flipud(t2.hs);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %       t2.tm = t2.tm(inInds3,:,:,:);t2.tm = flipud(t2.tm);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %       t2.tp = t2.tp(inInds3,:,:,:);t2.tp = flipud(t2.tp);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %       t2.x = t2.x(inInds3);t2.x = flipud(t2.x);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %       t2.y = t2.y(inInds3);t2.y = flipud(t2.y);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %       t2.x_u = t2.x_u(inInds3);t2.x_u = flipud(t2.x_u);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %       t2.y_u = t2.y_u(inInds3);t2.y_u = flipud(t2.y_u);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %       t2.lat = t2.lat(inInds3);t2.lat = flipud(t2.lat);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %       t2.lon = t2.lon(inInds3);t2.lon = flipud(t2.lon);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %       % Combine in to one structure
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %       L = struct;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %       L.hs = vertcat(t1.hs,t2.hs);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %       L.tm = vertcat(t1.tm,t2.tm);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %       L.tp = vertcat(t1.tp,t2.tp);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %       L.x = vertcat(t1.x,t2.x);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %       L.y = vertcat(t1.y,t2.y);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %       L.x_u = vertcat(t1.x_u,t2.x_u);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %       L.y_u = vertcat(t1.y_u,t2.y_u);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %       L.lon = vertcat(t1.lon,t2.lon);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %       L.lat = vertcat(t1.lat,t2.lat);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     else
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         L = L3;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         disp('Using LUT Zone 3')
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         tranInds = inInds3; % Subsample transect to be only part in basin
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         L.hs = L.hs(tranInds,:,:,:);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         L.tm = L.tm(tranInds,:,:,:);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         L.tp = L.tp(tranInds,:,:,:);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         L.x = L.x(tranInds);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         L.y = L.y(tranInds);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         L.x_u = L.x_u(tranInds);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         L.y_u = L.y_u(tranInds);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         L.lat = L.lat(tranInds);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         L.lon = L.lon(tranInds);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     end