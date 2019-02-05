%   Creates a hindcast for a specific transect
%   This script will create a hindcast along any specific transect created
%   in google earth.  The user determines what type of data they want to
%   use as a forcing.  The data can be obs, quantile corrected NNRP, or
%   NNRP.  This function should be called in a script and then the output
%   should be used to plot.

% NOTE THIS FUCKER TAKES 4 HOURS TO RUN
clearvars
addpath C:\Functions_Matlab\cbrewer

slr_val = 2.3;% The 50% 2100 scenario
ft2m = .3048;
slr = slr_val*ft2m; %[m]

% ----------------- Folder Locations of data ------------------------------
nnrp_fol = 'E:\Abbas\PS_COSMOS\Thesis_Modeling\Quantile_Correction\NNRP_WaterPointData\';
mask_fol = 'E:\Abbas\PS_COSMOS\Thesis_Modeling\KML\PugetSoundShorline\SanJuanShorelines\NNRPMasks\';
tran_fol = 'E:\Abbas\PS_COSMOS\Thesis_Modeling\KML\PugetSoundShorline\SanJuanShorelines\';
lut_fol = 'E:\Abbas\PS_COSMOS\Thesis_Modeling\LUT\Offshore_LUT_ExtractedTransects\SanJuan_IslandsLUT\';
tide_fol = 'E:\Abbas\Model_Met_Forcings\COOPS_tides\';

% ----------------- Key Value Pairs for Basin and NNRP --------------------
keys = {'LopezIsland_Mask1.kml','Orcas_Mask1.kml','SanJuan_Mask1.kml','SanJuan_Mask2.kml',...
    'SanJuan_Mask3.kml','Stuart_Mask1.kml','Waldron_Mask1.kml','Waldron_Mask2.kml'};
values = {};
% Sort the keys alphabetically 
[keys,I] = sort(keys);
values = values(I);

DM = containers.Map(keys,values);


% ------------------ Load LUT ---------------------------------------------
luts = dir([lut_fol '*.mat']);
for ii = 1:length(luts)
    L(ii) = load([lut_fol luts(ii).name]);
end

% ----------------- Load in single NNRP Point for the time vector ---------
nn2load = [nnrp_fol 'JDF1_NNRP_Point.mat'];
N = load(nn2load);
N.time  = N.time';
n_time = N.time;
% Interp NNRP from 6 hourly to hourly
time = datenum(n_time(1):1/24:n_time(end));
clear N

% ------------------ Load in Tides ----------------------------------------
tide_stns = {'tacoma\Tacoma_Reconstruct_NAVD88.mat';'seattle\seattle_hrV_NAVD88.mat'};
stn_nm = {'tacoma';'seattle'};
mhhw = [2.86,2.75];
for ii = 1:length(tide_stns)
    if ii == 1 % Tacoma has different format
        tides = load([tide_fol tide_stns{ii}]);
        T(ii).name = stn_nm{ii};
        T(ii).twl = tides.twl_patched;
        T(ii).time = tides.time;
    else
        load([tide_fol tide_stns{ii}]);
        T(ii).name = stn_nm{ii};
        T(ii).twl = tides.WL_VALUE;
        T(ii).time = tides.time;
    end
    % Interpolate twl onto time vec
    inds = find(unique(T(ii).time));
    T(ii).time = T(ii).time(inds);
    T(ii).twl = T(ii).twl(inds);
    inds = find(diff(T(ii).time) == 0);
    T(ii).time(inds) = [];
    T(ii).twl(inds) = [];
    T(ii).twl_i = double(interp1(T(ii).time,T(ii).twl,time));
    T(ii).MHHW = mhhw(ii);
end

% ----------------- Load Transect -----------------------------------------
% K = kml2struct([tran_fol 'Swan_10mContour.kml']);


disp('Done Loading Data - Beginning Hindcast')
%% Create Hindcast from LUT
% ------------- Dimensions of LUT -----------------------------------------
% ------------- setup for interp ------------------------------------------
[X,Y,Z] = meshgrid(L(1).tide,L(1).speed,L(1).direc);
n_inds = [];
% ------------- Loop over basins then extraction points -------------------
tic
for kk = 1:length(L) % Loop through the basin masks and grab subset of points to make hindcast of
    
    % ------------------- Load Forcings -----------------------------------
    % Grab NNRP Forcing
    NN_point = DM(keys{kk});
    N = load([nnrp_fol NN_point]);
    N.wnddir = squeeze(N.wnddir);
    N.wndspd = squeeze(N.wndspd);
    
    % Convert to U and V space
    [N.U,N.V] = wind2UV('from',N.wnddir,'compass',N.wndspd);
    U = interp1(N.time,N.U,time)';
    V = interp1(N.time,N.V,time)';
    % Convert Back to wnddir and speed
    speed = sqrt(U.^2+V.^2);
    wnddir = wrap2360(270-(atan2(V,U)*(180/pi)));
    
    % Grab Water level
    twl = M(kk).twl';
    
    % --------------- Downsample basin to be 100 point segments -----------
    % Due to space restraints with data, I will be saving 100 length
    % segments of each transect.  First find out how many segments you will
    % have.  Then loop through each segment and save.
    numFiles = floor(length(L.x)/100); % Number of complete 1000 point transects
    extraPoints = rem(length(L.x),100); % What is left over that isn't part of 1000
    inds = 1:1:99; % start with first 1000 points, kind of like a count variable
    
    
    % --------- Loop through the LUT and interp to make Hindcast ----------
    for ff = 1:numFiles+1
        
        % Initialize Variables
        hs_ts = zeros(length(inds),length(time));
        hs_ts_slr = hs_ts;
        tp_ts = hs_ts;
        tp_ts_slr = hs_ts;
        
        if ff == numFiles+1 % As long as we are not on the remiander iteration
            inds = prevInds; % Go back to previous inds because they've been wrongfully incremented 
            inds = inds(end)+1:1:inds(end)+1+extraPoints;
            hs_ts = zeros(length(inds),length(time));
            hs_ts_slr = hs_ts;
            tp_ts = hs_ts;
            tp_ts_slr = hs_ts;
        end
        % Loop through indices 
        for nn = 1:length(inds)
            temp = permute(squeeze(L.hs(inds(nn),:,:,:)),[2 1 3]);
            nan_inds = isnan(temp);
            temp(nan_inds) = 0;
            if ~isnan(temp(2,1))
                hs_ts(nn,:) = interp3(X,Y,Z,temp,twl,speed,wnddir,'linear')'; % Hsig timeseries
                hs_ts_slr(nn,:) = interp3(X,Y,Z,temp,twl+slr,speed,wnddir,'linear')';
                %---Peak Period
                temp = permute(squeeze(L.tp(inds(nn),:,:,:)),[2 1 3]);
                nan_inds = isnan(temp);
                temp(nan_inds) = 0;
                tp_ts(nn,:) = interp3(X,Y,Z,temp,twl,speed,wnddir,'linear')';
                tp_ts_slr(nn,:) = interp3(X,Y,Z,temp,twl+slr,speed,wnddir,'linear')';

                fprintf('Completed %d out of %d...Moving On\n',inds(nn),length(L.x))
            else
                fprintf('%d out of %d is a NAN point...Moving On',inds(nn),length(L.x))
            end
        end
        % ------------ Set any negative hs to zero (caused by winds < 2.5m/s) -----
        hs_ts(hs_ts<0) = 0;
        hs_ts_slr(hs_ts_slr<0) = 0;
        
        % ------------ Set any negative hs to zero (caused by winds < 2.5m/s) -----
        tp_ts(tp_ts<0) = 0;
        tp_ts_slr(tp_ts_slr<0) = 0;
        
        % Status Message
        fprintf('Done with %s: Zone %d of %d\n',basinName,ff,numFiles+1)
        
        % Save
        x_u = L.x_u(inds);
        y_u = L.y_u(inds);
        lon = L.lon(inds);
        lat = L.lat(inds);
        slrMeters = slr;
        slrFeet = slr_val;
        
        outName = sprintf('%s_%d_%d_waveOut.mat',basinName,ff,numFiles+1);
        save(outName,'x_u','y_u','lon','lat','hs_ts','hs_ts_slr','tp_ts',...
            'tp_ts_slr','speed','wnddir','twl','time','slrMeters','slrFeet')
        disp('Wave Hindcast Saved...Moving On')
        movefile(outName,'E:\Abbas\PS_COSMOS\Thesis_Modeling\LUT\Hindcast_Output')
        
        % Increment Inds to next group
        prevInds = inds;
        inds = inds(end)+1:1:inds(end)+100;
    end
    fprintf('%s Fully Completed - %.1f Percent of Basins Completed',basinName,round((kk/length(keys)*100),1))
end

toc