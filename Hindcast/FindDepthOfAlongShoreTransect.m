% Gets the depth at each alongshore transect point 
clearvars

% First load in the depths for the jdf and ps model 
% ------------------------- JDF Model -------------------------------------
dep_nm = 'C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\SWAN\PS_RegionalModel\INP_sph2\jdaf_new.BOT';
grd_nm = 'C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\SWAN\PS_RegionalModel\INP_sph2\jdf_sph_swn.grd';

% From .swn file
mx = 659;
my = 1195;
% Load grid
J = swan_io_grd('read',grd_nm,mx,my,99,99);

% Load bottom
S.mxinp = mx;
S.myinp = my;
S.idla = 4;
S.nhedf = 0;
S.fname1 = dep_nm;
S.quantity = 'depth';
J.Z = swan_io_bot('read',grd_nm,S);
% -------------------------- PS Model -------------------------------------

dep_nm = 'C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\SWAN\PS_RegionalModel\INP_sph\PugetSound2.BOT';
grd_nm = 'C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\SWAN\PS_RegionalModel\INP_sph\pug8_sph_swn.grd';

% From .swn file
mx = 1555;
my = 524;
% Load grid
P = swan_io_grd('read',grd_nm,mx,my,99,99);

% Load bottom
S.mxinp = mx;
S.myinp = my;
S.idla = 4;
S.nhedf = 0;
S.fname1 = dep_nm;
S.quantity = 'depth';
P.Z = swan_io_bot('read',grd_nm,S);

%% Load in Alongshore transect
kml_fol = 'C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\KML\PugetSoundShorline\'; % Location of shoreline KML to use
kml_mask = 'C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\KML\Basin_Masks\LargeBasinMasks\South_PS_SWAN_Domain.kml';
kml = 'Swan_10mContour.kml'; 

% Load KML Data - Basin Mask and Transect 
K = kml2struct([kml_fol kml]);  
M = kml2struct(kml_mask);

% Convert to utm
[K.x,K.y] = deg2utm(K.Lat,K.Lon);

% Create finely spaced alongshore transect to be derefined later 
x_e = [];
y_e = [];
for nn = 1:length(K.x)-1
step = 1; % 1 meter
[ line_x, line_y ] = createTransect( K.x(nn), K.y(nn), K.x(nn+1), K.y(nn+1), step );
x_e = cat(1,x_e,line_x');
y_e = cat(1,y_e,line_y');
end
% Get rid of any nans
nan_inds = isnan(x_e);
x_e(isnan(x_e)) = [];
y_e(isnan(y_e)) = [];

% Derefine to be specific resolution 
res = 50; %[m]
ra = 0:res:length(x_e);
ra(1) = 1;
x_e = x_e(ra);
y_e = y_e(ra);
% Convert back to lat lon to find grid indices 
[lat,lon] = utm2deg(x_e,y_e,char(ones(length(x_e),1)*'10 T'));

%% Now find points that are within the JDF, PS, and SGA domains 

% Find which points are within Southern Domain 
in = [];
inp = inpolygon(lon,lat,M.Lon,M.Lat); % Find what part of the transect is in the Basin
inds = find(inp == 1); % Find the basin points
P.lon = lon(inds);
P.lat = lat(inds);
P.x = x_e(inds);
P.y = y_e(inds);
in = vertcat(in,inds);

% Find which points are within JDF domain 
kml_mask = 'C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\KML\Basin_Masks\LargeBasinMasks\JDF_SWAN_Domain.kml'; % Used to determine which model to use
M = kml2struct(kml_mask);
inp = inpolygon(lon,lat,M.Lon,M.Lat); % Find what part of the transect is in the Basin
inds = find(inp == 1); % Find the basin points
J.lon = lon(inds);
J.lat = lat(inds);
J.x = x_e(inds);
J.y = y_e(inds);
in = vertcat(in,inds); % Add them to the list

% Remaining Inds must be the SGA domain 
inds = 1:1:in(1)-1;
G.lon = lon(inds);
G.lat = lat(inds);
G.x = x_e(inds);
G.y = y_e(inds);

%% Now I need to get the depths at each location of each 'Sub-domain'
for bb = 1:3
    if bb == 1
        domain = 'JDF';
    elseif bb == 2
        domain = 'PS';
    else
        domain = 'SGA';
    end
    
    switch domain
        case 'JDF'
            B = J;
        case 'PS'
            B = P;
        case 'SGA'
            B = J;
            B.lon = G.lon;
            B.lat = G.lat;
            B.x = G.x;
            B.y = G.y;
    end
    
    i_x = zeros(length(B.lon),1);
    i_y = i_x;
    % Find nearest grid cells for extraction
    for nn = 1:length(B.lon)
        [i_x(nn), i_y(nn)] = findNearestGridPoint(B.X,B.Y,B.lon(nn),B.lat(nn));
    end
    
    [n, m] = size(B.X);
    
    % Save output locations and grab the depths 
    B.xp = B.X(sub2ind([n m],i_x,i_y));
    B.yp = B.Y(sub2ind([n m],i_x,i_y));
    B.zp = B.Z(sub2ind([n m],i_x,i_y));
    switch domain
        case 'JDF'
            J.zp = B.zp;
            J.xp = B.xp;
            J.yp = B.yp;
        case 'PS'
            P.zp = B.zp;
            P.xp = B.zp;
            P.yp = B.yp;
        case 'SGA'
            G.zp = B.zp;
            G.xp = B.xp;
            G.yp = B.yp;
            % Make Along Transect 'Transect'
            s = sqrt(B.x.^2+B.x.^2);
            s = s - min(s);
            G.s = s;
    end
end

