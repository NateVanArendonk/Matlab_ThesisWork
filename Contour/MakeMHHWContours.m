% First load in the DEM
E = load('E:\Abbas\Modeling Resources\PS_DEM\Ruston_Way\RustonWayCONED_DEM.mat');
load('C:\Users\ahooshmand\Desktop\PS_COSMOS\Salish_Model_Resources\WA_Spatial_Data\WA_coast_UTM');

K = kml2struct(['OB_RW_Mask2.kml']);
[K.x,K.y] = deg2utm(K.Lat,K.Lon);
tic
IN = inpolygon(E.x,E.y,K.x,K.y);
toc
E.x = E.x(IN);
E.y = E.y(IN);
E.z = E.z(IN);

dxdy = 1;
mx = min(E.x):dxdy:max(E.x);
my = min(E.y):dxdy:max(E.y);
[DX,DY] = meshgrid(mx,my);
tic
DZ = griddata(E.x,E.y,E.z,DX,DY);
toc


%%
siteMHHW = 2.83; %[m] - NAVD
site2070 = 3.23; %[m] - NAVD - 50% scenario by 2070, 1.3 feet SLR
site2120 = 3.69; %[m] - NAVD - 50% scenario by 2120, 2.8 feet SLR


c = contourc(mx,my,DZ,[2.83 2.83]);
c(:,1) = []; 
x1 = c(1,:);
y1 = c(2,:);
bID = find(x1 < 1000);
x1 = x1(1:bID-1);
y1 = y1(1:bID-1);
%% Output to KML

utmzone = repmat('10 T',length(x1),1);
[lat,lon] = utm2deg(x1,y1,utmzone);

% Write KML
kmlwriteline('OBRW_MHHW.kml',lat,lon);