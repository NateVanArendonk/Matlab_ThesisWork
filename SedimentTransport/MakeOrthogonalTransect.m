clearvars
clc

% Load in Bathy/DEM
% File locations
dep_nm = 'E:\Abbas\PS_COSMOS\Thesis_Modeling\SWAN\PS_RegionalModel\INP_sph2\jdaf_new.BOT';
grd_nm = 'E:\Abbas\PS_COSMOS\Thesis_Modeling\SWAN\PS_RegionalModel\INP_sph2\jdf_sph_swn.grd';
% From .swn file
mxc = 659;
myc = 1195;
% Load grid
J = swan_io_grd('read',grd_nm,mxc,myc,99,99);

% Load bottom
S.mxinp = mxc;
S.myinp = myc;
S.idla = 4;
S.nhedf = 0;
S.fname1 = dep_nm;
S.quantity = 'depth';
J.Z = swan_io_bot('read',grd_nm,S);


% File locations
dep_nm = 'E:\Abbas\PS_COSMOS\Thesis_Modeling\SWAN\PS_RegionalModel\INP_sph\PugetSound2.BOT';
grd_nm = 'E:\Abbas\PS_COSMOS\Thesis_Modeling\SWAN\PS_RegionalModel\INP_sph\pug8_sph_swn.grd';
% From .swn file
mxc = 1555;
myc = 524;

% Load grid
P = swan_io_grd('read',grd_nm,mxc,myc,99,99);

% Load bottom
S.mxinp = mxc;
S.myinp = myc;
S.idla = 4;
S.nhedf = 0;
S.fname1 = dep_nm;
S.quantity = 'depth';
P.Z = swan_io_bot('read',grd_nm,S);
clear S mxc myc

%% Make a coarse bathy dem for sampling for direciton of shore
% Get DEM from both datasets
X = [J.X(:);P.X(:)];
Y = [J.Y(:);P.Y(:)];
Z = [J.Z(:);P.Z(:)];
% Get rid of NaN's
nInds = isnan(X);
X(nInds) = [];
Y(nInds) = [];
Z(nInds) = [];
% Convert to UTM
[X,Y] = deg2utm(Y,X);
% Make grid
dx = 250;
mx = min(X):dx:max(X);
my = min(Y):dx:max(Y);

% interp on to grid
[MX,MY] = meshgrid(mx,my);
MZ = griddata(X,Y,Z,MX,MY);

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
%%
clf
plot(x,y,'r')
hold on 
snX = cell(length(x)-2,1);
snY = snX;
snAz = zeros(length(x)-2,1);
tic
% Loop through and calculate perpendicular line
for ii = 2:length(x)-1
    % Get current point, and surrounding points for calculation of slope of normal line
    x_p = x(ii);
    y_p = y(ii);
    
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
    [line_x,line_y] = createShoreNormalOutward(m,x_p,y_p,MX,MY,MZ,50);
%     snX{ii} = line_x;
%     snY{ii} = line_y;
    
    plotting = 1;
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

for hh = 24
   W = load([hFol d(hh).name]);
   for ii = 1:length(W.lon)
       P = zeros(size(W.hs_ts));
       P2 = P;
       I = P;
       Sxy = P;
       
       % Find which location we are going to use 
       dist = sqrt((W.lon(ii) - lon).^2 + (W.lat(ii) - lat).^2);
       [~,Ii] = min(dist);
       if Ii == 1;
           Ii = 2;
       elseif Ii == length(lon)
           Ii = length(lon)-1;
       end
       az = snAz(Ii);
       line_x = snX{Ii};
       line_y = snY{Ii};
       
       % Subsample the Waves and Period  Direction 
       hs = W.hs_ts(ii,:);
       tm = W.tp_ts(ii,:); %%%%% WHEN YOU HAVE MEAN PERIOD SWAP IT OUT FOR THIS 
       wDir = W.dir_ts(ii,:); % Given in Nautical convention - Coming From 
       
       % Now calculate difference between waves and azimuth, normalizing to
       % a local coordinate system
       normalizeDeg = @(x)(-mod(-x+180,360)+180);
       DiffDeg = @(a,b)-1*(normalizeDeg(a) - normalizeDeg(b));
       thetaR = DiffDeg(az,wDir);
       
       % Calculate Sxy 
       % convert period to frequency 
       depth = 15; %meters, deep water (offshore) - arbitrary, assuming deep water waves for this domain 
       fo = 1./tm; % frequency
       % Calculate cp and cg 
       [cp, cg] = getcp(fo,depth); % cp is phase speed of wave, cg is group velocity 
       % Calculate Sxy
       Sxy(jj,:) = (cg./cp) .* ((sin(2*thetaR))/2);
       
       % Calculate I - Is this the alongshore current part?
       k = 0.4; % - Constant 
       I(jj,:) = k.*cp.*Sxy;
       
       % Estimate Wave Power 
       rho = 1024; % density of water 
       g = 9.81; % acceleration due to gravity 
       Ptemp = ((rho*(g^2))/64*pi).*(hs.^2).*tm;
       P(jj,:) = Ptemp./1000; % Convert to kiloWatts so divide by 1000
       
       P2(jj,:) = 0.5 .*(hs.^2).*tm;
   end
   
   
   
end
