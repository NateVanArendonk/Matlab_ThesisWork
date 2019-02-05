% Extract Wave Grids at Each Polygon
clearvars

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

% Load WCRP KMLs
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
end

kml_mask = 'E:\Abbas\PS_COSMOS\Thesis_Modeling\KML\Basin_Masks\LargeBasinMasks\JDF_SWAN_Domain.kml'; % Used to determine which model to use
J = kml2struct(kml_mask);
kml_mask = 'E:\Abbas\PS_COSMOS\Thesis_Modeling\KML\Basin_Masks\LargeBasinMasks\PS_SWAN_Domain.kml';
P = kml2struct(kml_mask);

%%  Pick which wave grid to choose from is the next step
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

% NNRP Points
inds = [];
np = dir('E:\Abbas\PS_COSMOS\Thesis_Modeling\Quantile_Correction\NNRP_PointData\*.mat');
for f = 1:length(np)
    f_open = strcat('E:\Abbas\PS_COSMOS\Thesis_Modeling\Quantile_Correction\NNRP_PointData\',np(f).name);
    n = load(f_open);
    %     periodInd = strfind(np(f).name);
    N(f).name = np(f).name;%(1:periodInd-1);
    N(f).lat = n.lat(1);
    N(f).lon = n.lon(1);
    [N(f).x,N(f).y] = deg2utm(n.lat(1),n.lon(1));
    N(f).land = n.land;
    if n.land
        inds(end+1) = f;
    end
end
% Get rid of Land Points
N(inds) = [];
for ii = 1:length(K)
    dist = zeros(length(N),1);
    for nn = 1:length(N)
        dist(nn) = sqrt((K(ii).cx - N(nn).x).^2 + (K(ii).cy - N(nn).y).^2);
    end
    [~,I] = min(dist);
    K(ii).NN_Point = I;
end

% Loop through and plot closest NNRP point to circle (NNRP Point is water
% point)
plotting = 0;
if plotting
    for mm = 1:length(K)
        clf
        for ii = 1:length(K)
            plot(K(ii).x,K(ii).y,'k')
            hold on
            plot(K(ii).cx,K(ii).cy,'r*')
        end
        axis equal
        for ii = 1:length(N)
            plot(N(ii).x,N(ii).y,'ob')
        end
        
        
        plot(K(mm).x,K(mm).y,'g')
        plot(K(mm).cx,K(mm).cy,'y')
        plot(N(K(mm).NN_Point).x,N(K(mm).NN_Point).y,'bo','MarkerFaceColor','b')
        pause(0.3)
        
    end
end

%% Loop through and find grid points within each polygon 
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