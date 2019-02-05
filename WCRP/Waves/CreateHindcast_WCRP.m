%   Creates a hindcast for a specific transect
%   This script will create a hindcast along any specific transect created
%   in google earth.  The user determines what type of data they want to
%   use as a forcing.  The data can be obs, quantile corrected NNRP, or
%   NNRP.  This function should be called in a script and then the output
%   should be used to plot.

clearvars


slr_val = 0;% The 50% 2100 scenario
ft2m = .3048;
slr = slr_val*ft2m; %[m]

% ----------------- Folder Locations of data ------------------------------
nnrp_fol = 'E:\Abbas\PS_COSMOS\Thesis_Modeling\Quantile_Correction\NNRP_PointData\';
lut_fol = 'E:\Abbas\WCRP\WaveGridExtraction\WCRP_LUTs\';
tide_fol = 'E:\Abbas\Model_Met_Forcings\COOPS_tides\';

% Get list of files in each Folder - NNRP and LUT 
Q = dir('E:\Abbas\WCRP\WaveGridExtraction\WCRP_LUTs\*.mat');
M = dir('E:\Abbas\PS_COSMOS\Thesis_Modeling\Quantile_Correction\NNRP_WaterPointData\*.mat');

% ----------------- Load in single NNRP Point for the time vector ---------
nn2load = [nnrp_fol 'Blaine_NNRP_Point.mat'];
N = load(nn2load);
N.time  = N.time';
n_time = N.time;
% Interp NNRP from 6 hourly to hourly
time = datenum(n_time(1):1/24:n_time(end));
clear N

% ------------------ Load in Tides ----------------------------------------
tide_stns = {'tacoma\Tacoma_Reconstruct_NAVD88.mat';'seattle\seattle_hrV_NAVD88.mat'};
stn_nm = {'tacoma';'seattle'};
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
end

% Load in a LUT to get dimensions for Interp
L1 = load([lut_fol, 'WCRP_Zone_-124.7_48.4.mat']);

% ----------------- Load in Masks -----------------------------------------
kml_mask = 'E:\Abbas\PS_COSMOS\Thesis_Modeling\KML\Basin_Masks\LargeBasinMasks\JDF_SWAN_Domain.kml'; % Used to determine which model to use
J = kml2struct(kml_mask);
kml_mask = 'E:\Abbas\PS_COSMOS\Thesis_Modeling\KML\Basin_Masks\LargeBasinMasks\PS_SWAN_Domain.kml';
P = kml2struct(kml_mask);

%% Load in NNRP Data

% NNRP Points
np = dir('E:\Abbas\PS_COSMOS\Thesis_Modeling\Quantile_Correction\NNRP_WaterPointData\*.mat');
for f = 1:length(np)
    f_open = strcat('E:\Abbas\PS_COSMOS\Thesis_Modeling\Quantile_Correction\NNRP_WaterPointData\',np(f).name);
    n = load(f_open);
    %     periodInd = strfind(np(f).name);
    N(f).name = np(f).name;%(1:periodInd-1);
    N(f).lat = n.lat(1);
    N(f).lon = n.lon(1);
    [N(f).x,N(f).y] = deg2utm(n.lat(1),n.lon(1));
end
disp('Done Loading Data - Beginning Hindcast')
%% Create Hindcast from LUT
% ------------- Dimensions of LUT -----------------------------------------
% ------------- setup for interp ------------------------------------------
[X,Y,Z] = meshgrid(L1.tide,L1.speed,L1.direc);
n_inds = [];
% ------------- Loop over basins then extraction points -------------------
tic
for kk = 1:length(Q) % Loop through the circles
    
    L = load([lut_fol, Q(kk).name]);
    matInd = strfind(Q(kk).name,'.mat');
    hindName = Q(kk).name(1:matInd-1);
    if ~exist(['WCRP_Hindcast\',hindName], 'dir')
        mkdir(['WCRP_Hindcast\',hindName])
    end
    
    % Find closest NNRP point to each wave grid cell in the circle
    nnList = zeros(length(L.x),1); % List to house all of the NNRP points to use 
    for kk = 1:length(L.x) % Go through each point in the Circle
        dist = zeros(length(N),1); % House the distance between said NNRP point and the current wave grid point 
        for nn = 1:length(N)
            dist(nn) = sqrt((L.x(kk) - N(nn).lon).^2 + (L.y(kk) - N(nn).lat).^2);
        end
        [~,I] = min(dist);
        nnList(kk) = I;
    end
    
    % Grab Water level - Just using seattles for now
    twl = T(2).twl_i';
    
    % --------------- Downsample basin to be 100 point segments -----------
    % Due to space restraints with data, I will be saving 100 length
    % segments of each transect.  First find out how many segments you will
    % have.  Then loop through each segment and save.
    numFiles = floor(length(L.x)/1000); % Number of complete 1000 point transects
    extraPoints = rem(length(L.x),1000); % What is left over that isn't part of 1000
    if numFiles < 1
        inds = 1:1:length(L.x);
    else
        inds = 1:1:999; % start with first 100 points, kind of like a count variable
    end

    % --------- Loop through the points in the Circle and interp to make Hindcast ----------
    for ff = 1:numFiles+1        
        % Initialize Variables
        hs_BAD = zeros(length(inds),length(time));
%         hs_ts_slr = hs_ts;
        tp_BAD = hs_BAD;
%         tp_ts_slr = hs_ts;
        
        if ff == numFiles+1 && ff ~= 1 % As long as we are not on the remiander iteration
            inds = prevInds; % Go back to previous inds because they've been wrongfully incremented 
            inds = inds(end)+1:1:inds(end)+1+extraPoints;
            hs_BAD = zeros(length(inds),length(time));
%             hs_ts_slr = hs_ts;
            tp_BAD = hs_BAD;
%             tp_ts_slr = hs_ts;
        end
        hs = zeros(size(hs_BAD));
        tp = zeros(size(tp_BAD));
        % Loop through indices 
        for nn = 1:length(inds)
            
            %%%%%%%%%%%%%%%%%%%%%% YOU ARE RIGHT HERE LOADING IN THE NNRP
            %%%%%%%%%%%%%%%%%%%%%% POINT FOR EACH GRID CELL IN THE
            %%%%%%%%%%%%%%%%%%%%%% MODEL!!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % ------------------- Load Forcings ---------------------------
            % Grab NNRP Forcing
            NN_point =
            N = load([nnrp_fol M(NN_point).name]);
            N.wnddir = squeeze(N.wnddir);
            N.wndspd = squeeze(N.wndspd);
            
            % Convert to U and V space
            [N.U,N.V] = wind2UV('from',N.wnddir,'compass',N.wndspd);
            U = interp1(N.time,N.U,time)';
            V = interp1(N.time,N.V,time)';
            % Convert Back to wnddir and speed
            speed = sqrt(U.^2+V.^2);
            wnddir = wrap2360(270-(atan2(V,U)*(180/pi)));
            
            
            temp = permute(squeeze(L.hs(inds(nn),:,:,:)),[2 1 3]);
            nan_inds = isnan(temp);
            temp(nan_inds) = 0;
            if ~isnan(temp(2,1))
                temp2 = interp3(X,Y,Z,temp,twl,speed,wnddir,'linear')';
                hs(nn,:) = temp2;
                %hs_tp(nn,1:length(hs_tp)) = temp2;
                %hs_tp(nn,:) = interp3(X,Y,Z,temp,twl,speed,wnddir,'linear')'; % Hsig timeseries
%                 hs_ts_slr(nn,:) = interp3(X,Y,Z,temp,twl+slr,speed,wnddir,'linear')';
                %---Peak Period
                temp = permute(squeeze(L.tp(inds(nn),:,:,:)),[2 1 3]);
                nan_inds = isnan(temp);
                temp(nan_inds) = 0;
                tp(nn,:) = interp3(X,Y,Z,temp,twl,speed,wnddir,'linear')';
%                 tp_ts_slr(nn,:) = interp3(X,Y,Z,temp,twl+slr,speed,wnddir,'linear')';

                fprintf('Completed %d out of %d...Moving On\n',inds(nn),length(L.x))
            else
                fprintf('%d out of %d is a NAN point...Moving On',inds(nn),length(L.x))
            end
        end
        % ------------ Set any negative hs to zero (caused by winds < 2.5m/s) -----
        hs(hs<0) = 0;
%         hs_ts_slr(hs_ts_slr<0) = 0;
        
        % ------------ Set any negative hs to zero (caused by winds < 2.5m/s) -----
        tp(tp<0) = 0;
%         tp_ts_slr(tp_ts_slr<0) = 0;
        
        % Status Message
        fprintf('Done with %d out of %d\n',ff,numFiles+1)
        
        % Save
        px = L.x(inds);
        py = L.y(inds);
        pz = L.z(inds);
        model = L.Model;
        NN = M(NN_point).name;
        IN = L.IN;
        centLat = L.centLat;
        centLon = L.centLon;
        centX = L.centX;
        centY = L.centY;
        circ_lat = L.circ_lat;
        circ_lon = L.circ_lon;
        circ_x = L.circ_x;
        circ_y = L.circ_y;
        depth = L.depth(inds);
        slrMeters = slr;
        slrFeet = slr_val;
        
        outName = sprintf('%s_%d_%d_waveOut.mat',hindName,ff,numFiles+1);
        save(outName,'px','py','pz','model','NN','IN','centLat','centLon',...
            'centX','centY','circ_lat','circ_lon','circ_x','circ_y','hs','tp',...
            'speed','wnddir','twl','time','slrMeters','slrFeet','depth','-v7.3')
        disp('Wave Hindcast Saved...Moving On')
        movefile(outName,['WCRP_Hindcast\',hindName])
        
        % Increment Inds to next group
        prevInds = inds;
        inds = inds(end)+1:1:inds(end)+1000;
    end
    fprintf('%s Fully Completed - %.1f Percent of Circles Completed',hindName,round((kk/length(Q)*100),1))
end

toc