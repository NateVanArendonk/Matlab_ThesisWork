% Develop LUT from Nested SWAN Runs
% - Code modified from Skagit LUT written by S.Crosby: N.V. 12/7/2017
% 1. Load  model run data
% 2. Determine grid cell or two from NNRP predictions to use
% 3. Create time series of wind speed, dir, and water level
% 4. Create interpolation from LUT given speed, dir and wl
% 5. Determine output points (or transects)
clearvars
addpath C:/Functions_Matlab


% SWAN runs location
fol_run = 'E:\Abbas\PS_COSMOS\Thesis_Modeling\SWAN\PS_RegionalModel\Model_Runs\NV_LUT\RES2\'; % Location of Swan Model Runs
kml_fol = 'C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\KML\PugetSoundShorline\'; % Location of shoreline KML to use
kml_mask = 'C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\KML\Basin_Masks\LargeBasinMasks\South_PS_SWAN_Domain.kml'; % Used to determine which model to use

% Parameter Ranges
tide = [-2,-1,0,1,2,3,5.5];
speed = 5:5:30;
direc = -10:10:370;

% set minimum depth
min_depth = .1; %[m];

% Load alongshore offshore transect
kml = 'Swan_10mContour.kml';
kml_load = [kml_fol kml]; % KML of 10m isobar
out_file = 'SWAN_10m_LUT_offshore_extract';

% Load KML Data
K = kml2struct(kml_load);
M = kml2struct(kml_mask);
clear kml_load

%%

% Create finely spaced alongshore transect to be derefined later
lon_e = [];
lat_e = [];
for nn = 1:length(K.Lon)-1
    step = .000014; % 1 meter,  dx = res_in_meters/73000, Note I used 91500 its half between 110000 and 73000
    [ line_lon, line_lat ] = createTransect( K.Lon(nn), K.Lat(nn), K.Lon(nn+1), K.Lat(nn+1), step );
    lon_e = cat(1,lon_e,line_lon');
    lat_e = cat(1,lat_e,line_lat');
end
% Get rid of any nans
nan_inds = isnan(lon_e);
lon_e(isnan(lon_e)) = [];
lat_e(isnan(lat_e)) = [];

% Derefine to be specific resolution
res = 50; %[m]
lon_s = lon_e(1:res:end);
lat_s = lat_e(1:res:end);

% Find which parts of kml transect are within the southern domain
[in,~] = inpolygon(lon_e,lat_e,M.Lon,M.Lat);
% Subsample kml to be just southern PS
lon_e = lon_e(in);
lat_e = lat_e(in);

% Initialize by loading a single run
tt = 1;
ss = 1;
dd = 1;
fname = sprintf('SpatialPS_s%02d_d%d_t%d.mat',speed(ss)*10,dir(dd),tide(tt)*10);
data = load([fol_run fname]);

% Find nearest grid cells for extraction
for nn = 1:length(lon_e)
    [i_x(nn), i_y(nn)] = findNearestGridPoint(data.Xp,data.Yp,lon_e(nn),lat_e(nn));
end

% Size of grid
[n, m] = size(data.Hsig);

% Save output locations
E.x = data.Xp(sub2ind([n m],i_x,i_y));
E.y = data.Yp(sub2ind([n m],i_x,i_y));

% Check and make sure that we grabbed a deep enough point, if not, grab a
% new point
% Check and make sure that we grabbed a deep enough point, if not, grab a
% new point
highEnd = 12;
lowEnd = 60;
for ii = 1:length(E.z)
    count = 0;
    % Original Values if this doesn't work below
    OGX = E.x(ii);
    OGY = E.y(ii);
    OGZ = E.z(ii);
    while E.z(ii) < highEnd || E.z(ii) > lowEnd
        depthGrid = zeros(3,3);
        xGrid = depthGrid;
        yGrid = depthGrid;
        ixGrid = depthGrid;
        iyGrid = depthGrid;
        % Order to work through wave grid
        xx = [-1 0 1 -1 0 1 -1 0 1];
        yy = [1 1 1 0 0 0 -1 -1 -1];
        % Grab depths at each point in box around current grid point
        % Note Depths are Positive!
        for jj = 1:3
            for kk = 1:length(xx)
                mxx = data.Xp(sub2ind([n, m],i_x(ii)+xx(kk),i_y(ii)+yy(kk)));
                myy = data.Yp(sub2ind([n, m],i_x(ii)+xx(kk),i_y(ii)+yy(kk)));
                mzz = data.Depth(sub2ind([n, m],i_x(ii)+xx(kk),i_y(ii)+yy(kk)));
                xGrid(2+xx(kk),2+yy(kk)) = mxx;
                yGrid(2+xx(kk),2+yy(kk)) = myy;
                depthGrid(2+xx(kk),2+yy(kk)) = mzz;
                ixGrid(2+xx(kk),2+yy(kk)) = i_x(ii)+xx(kk);
                iyGrid(2+xx(kk),2+yy(kk)) = i_y(ii)+yy(kk);
            end
        end
        % Grab the deepest depth point, and if it is less than 10, the loop
        % will stop, otherwise, it will keep on going
        % Lets first see if any are greater than or equal to the
        % criteria I set
        logMat = depthGrid >= highEnd;
        if sum(logMat(:)) >= 1
            depthGrid(~logMat) = 999;
            diffMat = ones(size(depthGrid))*highEnd;
            diffMat = abs(depthGrid - diffMat);
            [~,Id] = min(diffMat(:),[],'omitnan');
            % Get Row and Column
            [I_row, I_col] = ind2sub(size(depthGrid),Id);
        else % we need to relocate - find the closest point to my criteria
            diffMat = ones(size(depthGrid))*highEnd;
            diffMat = abs(depthGrid - diffMat);
            [~,Id] = min(diffMat(:),[],'omitnan');
            % Get Row and Column
            [I_row, I_col] = ind2sub(size(depthGrid),Id);
        end
        % Change wave grid point to be new depth location
        E.x(ii) = xGrid(I_row,I_col);
        E.y(ii) = yGrid(I_row,I_col);
        E.z(ii) = depthGrid(I_row,I_col);
        i_x(ii) = ixGrid(I_row,I_col);
        i_y(ii) = iyGrid(I_row,I_col);
        count = count + 1;
        if count > 20
            %                 if E.z(ii) > OGZ
            %                     E.x(ii) = OGX;
            %                     E.y(ii) = OGY;
            %                     E.z(ii) = OGZ;
            %                 end
            if E.z(ii) > lowEnd
                fprintf('WARNING - DEPTH IS GREATER THAN USER THRESHOLD')
            end
            break
        end
    end
end

% Loop over parameterizations (takes 15 minutes)
tic
E.hs = NaN(length(i_x),length(tide),length(speed),length(dir));
E.tp = E.hs;
E.tm = E.hs;


for tt = 1:length(tide)
    for ss = 1:length(speed)
        for dd = 1:length(direc)
            
            fname = sprintf('SpatialPS_s%02d_d%d_t%d.mat',speed(ss)*10,dir(dd),tide(tt)*10);
            if direc(dd) == -10
                ind = find(direc == 350);
                fname = sprintf('SpatialPS_s%02d_d%d_t%d.mat',speed(ss)*10,direc(ind),tide(tt)*10);
            elseif direc(dd) == 360
                ind = find(direc == 0);
                fname = sprintf('SpatialPS_s%02d_d%d_t%d.mat',speed(ss)*10,direc(ind),tide(tt)*10);
            elseif direc(dd) == 370
                ind = find(direc == 10);
                fname = sprintf('SpatialPS_s%02d_d%d_t%d.mat',speed(ss)*10,direc(ind),tide(tt)*10);
            end
            
            
            data = load([fol_run fname]);
            dep = data.Depth(sub2ind([n m],i_x,i_y))';
            dry_inds = dep < min_depth;
            
            % Get alongshore point wave info
            E.hs(:,tt,ss,dd) = data.Hsig(sub2ind([n m],i_x,i_y))';
            E.tp(:,tt,ss,dd) = data.RTpeak(sub2ind([n m],i_x,i_y))';
            E.tm(:,tt,ss,dd) = data.Tm01(sub2ind([n m],i_x,i_y))';
            
            % Set NaN dry indices
            E.hs(dry_inds,tt,ss,dd) = 0;
            E.tp(dry_inds,tt,ss,dd) = 0;
            E.tm(dry_inds,tt,ss,dd) = 0;
            
        end
    end
    fprintf('Completed Tide level %d - Moving On\n',tide(tt))
end
toc

E.tide = tide;
E.speed = speed;
E.dir = dir;
E.tranX = lon_e;
E.tranY = lat_e;
E.gridX = data.Xp;
E.gridY = data.Yp;

save(out_file,'-struct','E','-v7.3')




