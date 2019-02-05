% Extract transect from asc file
% Extracts all transects in kml and saves as matlab files
clearvars

% Load bathy (takes 2-3 minutes)
tic
B = ascii2xyz('nkmos1.asc');
toc

% Load kml
T = kml2struct('XB_Transects.kml');

% Conver to UTM
for ii=1:length(T)
    [T(ii).x,T(ii).y] = deg2utm(T(ii).Lat,T(ii).Lon);
end

% Select Transect
for ii = 1:length(T)

% Set resolution of transect
ds = 1; %[m]

% Create transect
[ line_x, line_y ] = createTransect( T(ii).x(1), T(ii).y(1), T(ii).x(2), T(ii).y(2), ds );

% Extract subset of depths that covers transect
inds = B(:,1) >= min(line_x) & B(:,1) <= max(line_x) & B(:,2) >= min(line_y) & B(:,2) <= max(line_y);
bathy_x = B(inds,1);
bathy_y = B(inds,2);
bathy_z = B(inds,3);

% Find nearest depth for each location on transect 
tic
line_z = zeros(size(line_x));
for jj = 1:length(line_x)
    dist2 = (line_x(jj)-bathy_x).^2 + (line_y(jj)-bathy_y).^2;
    [mindist, I] = min(dist2);
    line_z(jj) = bathy_z(I);
    %disp(sqrt(mindist))
end
toc

% Reverse transect if needed, want ocean to be first
if line_z(1) > line_z(2)
   line_x = fliplr(line_x);
   line_y = fliplr(line_y);
   line_z = fliplr(line_z);
end

% Include alongtransect
s = 0:ds:ds*(length(line_x)-1);

% Save transect
fname = T(ii).Name;
L.x = line_x;
L.y = line_y;
L.z = line_z;
L.ds = ds;
L.s = s;
save([fname '_xyz.mat'],'-struct','L')

end













