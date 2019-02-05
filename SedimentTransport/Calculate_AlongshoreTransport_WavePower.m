clearvars
clc

%% First Load in SWAN grids and bathy

% Load in Bathy/DEM
% File locations
% dep_nm = 'E:\Abbas\PS_COSMOS\Thesis_Modeling\SWAN\PS_RegionalModel\INP_sph2\jdaf_new.BOT';
% grd_nm = 'E:\Abbas\PS_COSMOS\Thesis_Modeling\SWAN\PS_RegionalModel\INP_sph2\jdf_sph_swn.grd';
% % From .swn file
% mxc = 659;
% myc = 1195;
% % Load grid
% J = swan_io_grd('read',grd_nm,mxc,myc,99,99);
%
% % Load bottom
% S.mxinp = mxc;
% S.myinp = myc;
% S.idla = 4;
% S.nhedf = 0;
% S.fname1 = dep_nm;
% S.quantity = 'depth';
% J.Z = swan_io_bot('read',grd_nm,S);
%
%
% % File locations
% dep_nm = 'E:\Abbas\PS_COSMOS\Thesis_Modeling\SWAN\PS_RegionalModel\INP_sph\PugetSound2.BOT';
% grd_nm = 'E:\Abbas\PS_COSMOS\Thesis_Modeling\SWAN\PS_RegionalModel\INP_sph\pug8_sph_swn.grd';
% % From .swn file
% mxc = 1555;
% myc = 524;
%
% % Load grid
% P = swan_io_grd('read',grd_nm,mxc,myc,99,99);
%
% % Load bottom
% S.mxinp = mxc;
% S.myinp = myc;
% S.idla = 4;
% S.nhedf = 0;
% S.fname1 = dep_nm;
% S.quantity = 'depth';
% P.Z = swan_io_bot('read',grd_nm,S);
% clear S mxc myc

%% Make a coarse bathy dem for sampling for direciton of shore normal
% Get DEM from both datasets
% X = [J.X(:);P.X(:)];
% Y = [J.Y(:);P.Y(:)];
% Z = [J.Z(:);P.Z(:)];
% % Get rid of NaN's
% nInds = isnan(X);
% X(nInds) = [];
% Y(nInds) = [];
% Z(nInds) = [];
% % Convert to UTM
% [X,Y] = deg2utm(Y,X);
% % Make grid
% dx = 250;
% mx = min(X):dx:max(X);
% my = min(Y):dx:max(Y);
%
% % interp on to grid
% [MX,MY] = meshgrid(mx,my);
% MZ = griddata(X,Y,Z,MX,MY);

%% Load alongshore transect from LUT
% reorder lut to be in correct order as well
lutFol = 'E:\Abbas\PS_COSMOS\Thesis_Modeling\LUT\Offshore_LUT_Extracted\';
fileNames = {'SWAN_10m_JDFLUT_offshore_extract.mat';'SWAN_10m_PSLUT_offshore_extract.mat';'SWAN_10m_SGALUT_offshore_extract.mat'};
lat = [];
lon = [];
for ii = 1:length(fileNames)
    temp = load([lutFol fileNames{ii}]);
    tlat = flipud(temp.lat);
    tlon = flipud(temp.lon);
    lat = vertcat(lat,tlat);
    lon = vertcat(lon,tlon);
end
clear tlon tlat temp lutFol grd_nm fileNames dep_nm ii

% -----------------------------------------------------------------
% Convert from KML to UTM
[x,y] = deg2utm(lat,lon);
%% Calculate shore normal direction


snX = cell(length(x)-2,1);
snY = snX;
snAz = zeros(length(x)-2,1);
tic
% Loop through and calculate perpendicular line
for ii = 2:length(x)-1
    % Get current point, and surrounding points for calculation of slope of normal line
    x_p = x(ii);
    y_p = y(ii);
    
    % Neighboring Points
    x_neigh1 = x(ii-1);
    y_neigh1 = y(ii-1);
    [ny1,nx1] = utm2deg(x_neigh1,y_neigh1,'10 T');
    x_neigh2 = x(ii+1);
    y_neigh2 = y(ii+1);
    [ny2,nx2] = utm2deg(x_neigh2,y_neigh2,'10 T');
    
    % Calculate heading of line made between neightbors
    [~,heading] = distance(ny1,nx1,ny2,nx2);
    
    % Loop through and grab average slope between first point and
    % subsequent point in line between transects
    slope = zeros(2,1);
    slope(1,1) = (y(ii-1) - y(ii))/(x(ii-1) - x(ii)); % slope of previous point with current point
    slope(2,1) = (y(ii+1) - y(ii))/(x(ii+1) - x(ii)); % Slope of next point with current point
    
    % Get rid of any nans or zeros
    slope(isnan(slope)) = [];
    slope(slope == 0) = [];
    
    % Get average slope of all points
    slope = mean(slope);
    
    % Calculate slope of normal line
    m = -1/slope;
    
    % Create orthogonal transect
    [line_x,line_y] = createShoreNormalOutwardFAST(m,x_p,y_p,heading,50);
    
    plotting = 0;
    if plotting
        plot(line_x,line_y,'k')
    end
    
    % Calculating heading of normal line now
    utmzone = repmat('10 T',length(line_x),1);
    [Lat,Lon] = utm2deg(line_x,line_y,utmzone);
    snX{ii} = Lon;
    snY{ii} = Lat;
    [~,az] = distance(Lat(1),Lon(1),Lat(end),Lon(end));
    snAz(ii) = az;
end
toc

%% Load in the waves
hFol = 'E:\Abbas\PS_COSMOS\Thesis_Modeling\LUT\Hindcast_Output\';
d = dir([hFol '*.mat']);

basins2do = {'Anacortes','BellinghamBay','Blaine','Hood_Hook',...
    'Hood_Kitsap','Middle_Hood','Port_Ludlow','Nor_Hood',...
    'LummiBay','NorthCentral','Olympia','OlympiaNorth','OlympiaWest',...
    'PadillaBay','PortTownsend_Hood','Skagit','SouthCentral',...
    'Tacoma','Whidbey_JDF','WhidbeyAnacortes','WhidbeyBasin'...
    'JDF1','JDF2','JDF3','JDF4','JDF5','JDF6','JDF7',...
    'JDF8','JDF9','JDF10','JDF11'...
    'BoundaryBay_Mask','PointRoberts_Mask','Tsawwassen_Mask',...
    'RichmondBC_Mask','IonaBeachBC_Mask','BurrardInlet_Mask'};
tic
for hh = 1:length(d)
    
    % Get just the basin name
    bName = zeros(length(d(hh).name),1);
    for bb = 1:length(bName)
        curVal = d(hh).name(bb);
        if strcmp(curVal,'i') % will make it imaginary and muck it all up
            curVal = 'A';
        end
        bName(bb) = str2double(curVal);
    end
    n_inds = find(~isnan(bName));
    basinName = d(hh).name(1:n_inds(1)-2);
    
    if ismember(basinName,basins2do)
        ind = strfind(d(hh).name,'waveOut.mat');
        saveNm = d(hh).name(1:ind-1);
        saveNm = [saveNm 'alongshoreOut.mat'];
        if ~exist(['Alongshore_Out/' saveNm])
            W = load([hFol d(hh).name]);
            W = rmfield(W,'wnddir');
            W = rmfield(W,'tp_ts');
            W = rmfield(W,'time');
            W = rmfield(W,'speed');
            
            P = zeros(size(W.hs_ts));
            %             P2 = P;
            I = P;
            Sxy = P;
            
            % Find which location we are going to use for the azimuth and lines
            dist = sqrt((W.lon - lon).^2 + (W.lat - lat).^2);
            [~,Ii] = min(dist);
            if Ii == 1
                Ii = 2;
            elseif Ii == length(lon)
                Ii = length(lon)-1;
            end
            az = snAz(Ii);
            line_x = snX{Ii};
            line_y = snY{Ii};
            
            %         % Subsample the Waves and Period  Direction
            hs = W.hs_ts;
            tm = W.tm_ts;
            wDir = W.dir_ts; % Given in Nautical convention - Coming From
            
            % Now calculate difference between waves and azimuth, normalizing to
            % a local coordinate system
            thetaR = zeros(1,length(wDir));
            for aa = 1:length(wDir)
                thetaR(aa) = smallestSignedAngleBetween(az,wDir(aa));
            end
            
            % Calculate Sxy
            % convert period to frequency
            depth = 15; %meters, deep water (offshore) - arbitrary, assuming deep water waves for this domain
            fo = 1./tm; % frequency
            % Calculate cp and cg
            [cp, cg] = getcp(fo,depth); % cp is phase speed of wave, cg is group velocity
            % Calculate Sxy
            Sxy = (hs.^2) .* (cg./cp) .* ((sind(2.*thetaR))/2);
            
            % Get rid of times where the wind is offshore
            g1 = thetaR <= -95;
            g2 = thetaR >= 95;
            Sxy(g1) = 0;
            Sxy(g2) = 0;
            
            % Calculate I - Is this the alongshore component
            k = 0.4; % - Constant
            I = k.*cp.*Sxy;
            
            % Estimate Wave Power
            rho = 1024; % density of water
            g = 9.81; % acceleration due to gravity
            Ptemp = ((rho*g^2)/(64*pi)) .* (hs.^2) .* tm;
            P = Ptemp./1000; % Convert to kiloWatts so divide by 1000
            P(g1) = 0;
            P(g2) = 0;
            % Other formula
            %             P2(ii,:) = 0.5 .*(hs.^2).*tm;
            clear fo Ptemp cp cg
            
            
            Lat = W.lat;
            Lon = W.lon;
            save(saveNm,'P','I','Sxy','hs','tm','wDir','Lon','Lat','az','line_x','line_y','thetaR')
            movefile(saveNm,'Alongshore_Out')
            clear W P I Sxy wDir hs tm Lat Lon
            fprintf('Completed %s\n',saveNm)
        end
    end
    if rem(hh,100) == 0
        fprintf('%2.1f Percent Complete\n',100*(hh/(length(d))))
    end
    
end
toc
return
%%
[hs_sorted,IS] = sort(hs,'descend');
nanInds = isnan(hs_sorted);
hs_sorted(nanInds) = [];
IS(nanInds) = [];
maxDir = wDir(IS);
clf
plot(hs_sorted,maxDir,'.')
hold on
line([0 3],[wrapTo360(az - 90) wrapTo360(az - 90)],'Color','r')
line([0 3],[wrapTo360(az + 90) wrapTo360(az + 90)],'Color','r')
line([0 3],[az az],'Color','k','linestyle','--')
xlabel('Wave Height [m]')
ylabel('Wind Direction [degrees]')
set(gca,'FontSize',14)

a1 = fill([0 0 3 3],[0 wrapTo360(az+90) wrapTo360(az+90) 0],[.7 .7 .7],'LineStyle','none');
a1.FaceAlpha = 0.3;
a2 = fill([0 0 3 3],[wrapTo360(az-90) 360 360 wrapTo360(az-90)],[.7 .7 .7],'LineStyle','none');
a2.FaceAlpha = 0.3;
grid on
printFig(gcf,'Tacoma_Example_WavesAndDirection',[11 8.5],'png',300)