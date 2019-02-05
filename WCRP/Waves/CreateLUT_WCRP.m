% Develop LUT from Nested SWAN Runs for each circle WCRP
addpath C:/Functions_Matlab


% SWAN runs location
jdf_run = 'E:\Abbas\PS_COSMOS\Thesis_Modeling\SWAN\PS_RegionalModel\Model_Runs\NV_LUT\RES1\'; % Location of Swan Model Runs
ps_run = 'E:\Abbas\PS_COSMOS\Thesis_Modeling\SWAN\PS_RegionalModel\Model_Runs\NV_LUT\RES2\';

% Parameter Ranges
tide = [-2,-1,0,1,2,3,4,5.5];
speed = 0:5:30;
direc = -10:10:370;

% set minimum depth
min_depth = .1; %[m];

% First load in the depths for the jdf and ps model
% ------------------------- JDF Model -------------------------------------
dep_nm = 'E:\Abbas\PS_COSMOS\Thesis_Modeling\SWAN\PS_RegionalModel\INP_sph2\jdaf_new.BOT';
grd_nm = 'E:\Abbas\PS_COSMOS\Thesis_Modeling\SWAN\PS_RegionalModel\INP_sph2\jdf_sph_swn.grd';

% From .swn file
mx = 659;
my = 1195;
% Load grid
JDF = swan_io_grd('read',grd_nm,mx,my,99,99);

% Load bottom
S.mxinp = mx;
S.myinp = my;
S.idla = 4;
S.nhedf = 0;
S.fname1 = dep_nm;
S.quantity = 'depth';
JDF.Z = swan_io_bot('read',grd_nm,S);

% -------------------------- PS Model -------------------------------------
dep_nm = 'E:\Abbas\PS_COSMOS\Thesis_Modeling\SWAN\PS_RegionalModel\INP_sph\PugetSound2.BOT';
grd_nm = 'E:\Abbas\PS_COSMOS\Thesis_Modeling\SWAN\PS_RegionalModel\INP_sph\pug8_sph_swn.grd';

% From .swn file
mx = 1555;
my = 524;
% Load grid
PS = swan_io_grd('read',grd_nm,mx,my,99,99);

% Load bottom
S.mxinp = mx;
S.myinp = my;
S.idla = 4;
S.nhedf = 0;
S.fname1 = dep_nm;
S.quantity = 'depth';
PS.Z = swan_io_bot('read',grd_nm,S);

% --- Load WCRP KMLs
kml_fol = 'E:\Abbas\PS_COSMOS\Thesis_Modeling\KML\WCRP_KML\';
kml_name = 'wcrp_locations.kml';

temp = kml2struct([kml_fol kml_name]);
C = kml2struct([kml_fol 'Circles2Cut.kml']);


for ii = 1:length(temp)
    if ~inpolygon(temp(ii).Lon,temp(ii).Lat,C.Lon,C.Lat)
        %         pgon = polyshape([temp(ii).Lon],[temp(ii).Lat]);
        %         plot(pgon);
        %         hold on
        K(ii).lon = temp(ii).Lon;
        K(ii).lat = temp(ii).Lat;
        K(ii).lon(end) = [];
        K(ii).lat(end) = [];
        [K(ii).x,K(ii).y] = deg2utm(K(ii).lat,K(ii).lon);
    else
        K(ii).lat = NaN;
        K(ii).lon = NaN;
        K(ii).x = NaN;
        K(ii).y = NaN;
    end
end
inds = [];
for ii = 1:length(K)
    if isnan(K(ii).lat)
        inds(end+1) = ii;
    end
end
% Get rid of shapes not part of my project
K(inds) = [];

% Add Center of circle to structure
for ii = 1:length(K)
    [K(ii).cx,K(ii).cy,~] = CircleFitByPratt([K(ii).x,K(ii).y]);
    [K(ii).ly,K(ii).lx] = utm2deg(K(ii).cx,K(ii).cy,'10 T');
end

kml_mask = 'E:\Abbas\PS_COSMOS\Thesis_Modeling\KML\Basin_Masks\LargeBasinMasks\JDF_SWAN_Domain.kml'; % Used to determine which model to use
J = kml2struct(kml_mask);
kml_mask = 'E:\Abbas\PS_COSMOS\Thesis_Modeling\KML\Basin_Masks\LargeBasinMasks\PS_SWAN_Domain.kml';
P = kml2struct(kml_mask); 
return
%%  Pick which wave grid to choose from 
Jinds = [];
Pinds = [];
Ginds = [];
for ii = 1:length(K)
    inds = inpolygon(K(ii).lon,K(ii).lat,J.Lon,J.Lat);
    tot = sum(inds);
    if tot >= length(K(ii).lon)/2 % If more than half of poly is within JDF
        Jinds(end+1) = ii; % add it to JDF tally
    else % Otherwise see how much is within PS model
        inds = inpolygon(K(ii).lon,K(ii).lat,P.Lon,P.Lat);
        tot = sum(inds);
        if tot >= length(K(ii).lon)/2% If more than half is within PS
            Pinds(end+1) = ii; % add it to PS tally
        else % otherwise it'll be in the Georgia Model
            Ginds(end+1) = ii;
        end
    end
end
% Combine Indices for JDF and Georgia
Jinds = [Jinds Ginds];

%% Load in NNRP Data

% % NNRP Points
% inds = [];
% np = dir('E:\Abbas\PS_COSMOS\Thesis_Modeling\Quantile_Correction\NNRP_WaterPointData\*.mat');
% for f = 1:length(np)
%     f_open = strcat('E:\Abbas\PS_COSMOS\Thesis_Modeling\Quantile_Correction\NNRP_WaterPointData\',np(f).name);
%     n = load(f_open);
%     %     periodInd = strfind(np(f).name);
%     N(f).name = np(f).name;%(1:periodInd-1);
%     N(f).lat = n.lat(1);
%     N(f).lon = n.lon(1);
%     [N(f).x,N(f).y] = deg2utm(n.lat(1),n.lon(1));
% %     N(f).land = n.land;
% %     if n.land
% %         inds(end+1) = f;
% %     end
% end
% % % Get rid of Land Points
% % N(inds) = [];
% % % Find closest wet NNRP point to the circle center 
% % for ii = 1:length(K)
% %     dist = zeros(length(N),1);
% %     for nn = 1:length(N)
% %         dist(nn) = sqrt((K(ii).cx - N(nn).x).^2 + (K(ii).cy - N(nn).y).^2);
% %     end
% %     [~,I] = min(dist);
% %     K(ii).NN_Point = I;
% % end

%% Loop through and find grid points within each circle polygon and store those indices 
for ii = 1:length(K)
    in = ismember(ii,Jinds);
    if in
        K(ii).Model = 'JDF';
        IN = inpolygon(JDF.X,JDF.Y,K(ii).lon,K(ii).lat);
    else
        K(ii).Model = 'PS';
        IN = inpolygon(PS.X,PS.Y,K(ii).lon,K(ii).lat);
    end
    K(ii).IN = {IN};
end

%% Build LUT

for ii = 1:length(K)
    tt = 1;
    ss = 3;
    dd = 10;
    tic
    
%     % Find closest NNRP point to each wave grid cell in the circle 
%     nnList = zeros(length(K(ii).lat),1);
%     for kk = 1:length(K(ii).lat)
%         dist = zeros(length(nnList),1);
%         for nn = 1:length(N)
%             dist(nn) = sqrt((K(ii).lon(kk) - N(nn).lon).^2 + (K(ii).lat(kk) - N(nn).lat).^2);
%         end
%         [~,I] = min(dist);
%         nnList(kk) = I;
%     end
    
    % Check and see if file already exists, if it does, move on to next LUT
    % circle 
    out_file = sprintf('WCRP_Zone_%3.1f_%2.1f.mat',K(ii).lx,K(ii).ly);
    if ~exist(['WCRP_LUTs\' out_file], 'file')
        
        if strcmp(K(ii).Model,'JDF') % If we are to use JDF Load a JDF Run
            fname = sprintf('SpatialJDF_s%02d_d%d_t%d.mat',speed(ss)*10,direc(dd),tide(tt)*10);
            fstart = 'SpatialJDF';
            fol_run = jdf_run;
        else % Otherwise Load a PS Run
            fname = sprintf('SpatialPS_s%02d_d%d_t%d.mat',speed(ss)*10,direc(dd),tide(tt)*10);
            fstart = 'SpatialPS';
            fol_run = ps_run;
        end
        data = load([fol_run fname]); % Load that Run
        
        % Get Depths Specific Grid of interest
        if strcmp(K(ii).Model,'JDF')
            data.Zp = JDF.Z';
        else
            data.Zp = PS.Z';
        end
        
        % reorient Grid to be same as inpolygon grid
        data.Xp = data.Xp';
        data.Yp = data.Yp';
        data.Botlev = data.Botlev';
        data.Depth = data.Depth';
        data.Hsig = data.Hsig';
        data.RTpeak = data.RTpeak';
        data.Tm01 = data.Tm01';
        
        % Size of grid
        [n, m] = size(data.Hsig);
        
        % Save output locations
        IN = K(ii).IN{1};
        E.x = data.Xp(IN);
        E.y = data.Yp(IN);
        E.z = data.Botlev(IN);
          
        % Loop over parameterizations (takes 15 minutes)
        tic
        E.hs = NaN(length(E.x),length(tide),length(speed),length(dir));
        E.tp = E.hs;
        E.tm = E.hs;
        
        for tt = 1:length(tide)
            for ss = 1:length(speed)
                for dd = 1:length(direc)
                    if speed(ss) ~= 0
                        fname = sprintf('%s_s%02d_d%d_t%d.mat',fstart,speed(ss)*10,direc(dd),tide(tt)*10);
                        if direc(dd) == -10
                            ind = find(direc == 350);
                            fname = sprintf('%s_s%02d_d%d_t%d.mat',fstart,speed(ss)*10,direc(ind),tide(tt)*10);
                        elseif direc(dd) == 360
                            ind = find(direc == 0);
                            fname = sprintf('%s_s%02d_d%d_t%d.mat',fstart,speed(ss)*10,direc(ind),tide(tt)*10);
                        elseif direc(dd) == 370
                            ind = find(direc == 10);
                            fname = sprintf('%s_s%02d_d%d_t%d.mat',fstart,speed(ss)*10,direc(ind),tide(tt)*10);
                        end
                        data = load([fol_run fname]);
                        % reorient Grid to be same as inpolygon grid
                        data.Xp = data.Xp';
                        data.Yp = data.Yp';
                        data.Depth = data.Depth';
                        data.Botlev = data.Botlev';
                        data.Hsig = data.Hsig';
                        data.RTpeak = data.RTpeak';
                        data.Tm01 = data.Tm01';
                        dry_inds = data.Depth(IN) < min_depth;  % Reverse depth, positive values = deeeeeeep
                        
                        % Get alongshore point wave info
                        E.hs(:,tt,ss,dd) = data.Hsig(IN);
                        E.tp(:,tt,ss,dd) = data.RTpeak(IN);
                        E.tm(:,tt,ss,dd) = data.Tm01(IN);
                        
                        % Set NaN dry indices
                        E.hs(dry_inds,tt,ss,dd) = 0;
                        E.tp(dry_inds,tt,ss,dd) = 0;
                        E.tm(dry_inds,tt,ss,dd) = 0;
                    else % Set values to zero if wind speed is less than 1
                        % Get alongshore point wave info
                        E.hs(:,tt,ss,dd) = 0;
                        E.tp(:,tt,ss,dd) = 0;
                        E.tm(:,tt,ss,dd) = 0;
                    end
                end
            end
            fprintf('Completed Tide level %d - Moving On\n',tide(tt))
        end
        toc
        
        % Save
        E.tide = tide; % WL
        E.speed = speed; % Windspeed
        E.direc = direc; % Wind Direction 
        E.Model = K(ii).Model; % Model that was used 
        E.IN = K(ii).IN; % Points in model that were used 
        E.NN_Point = K(ii).NN_Point; % NNRP point to use 
        E.circ_lat = K(ii).lat; % Latitude of WCRP Circle
        E.circ_lon = K(ii).lon; % Longitude of WCRP Circle 
        E.circ_x = K(ii).x; % UTM 'X' of WCRP Circle 
        E.circ_y = K(ii).y; % UTM 'Y' of WCRP Circle 
        E.centLat = K(ii).ly; % Latitude of Circle Center 
        E.centLon = K(ii).lx; % Longitude of Circle Center 
        E.centY = K(ii).cy; % UTM 'Y' of Circle Center 
        E.centX = K(ii).cx; % UTM 'X' of Circle Center 
        E.depth = data.Depth(IN);
        
        % Wave height and Peak/Mean period are already part of structure E
        
        save(out_file,'-struct','E','-v7.3')
        movefile(out_file,'WCRP_LUTs')
    end
    
    fprintf('----------------------------------------------------------------\n')
    fprintf('Completed Zone %d out of %d\n',ii,length(K))
    fprintf('----------------------------------------------------------------\n')
end
