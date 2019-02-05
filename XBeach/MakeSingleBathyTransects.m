clearvars
addpath C:\Functions_Matlab
addpath C:\Functions_Matlab\mapping\kml

% ---------- Folder of DEM
dem_fol = 'E:\Abbas\Modeling Resources\PS_DEM\Ruston_Way\CONED\';
fname = 'ruston_ascii.asc';

% ---------- Limits for Subsetting DEM
lim.lx = 5.33*10^5;
lim.ly = 5.2385*10^6;
lim.rx = 5.375*10^5;
lim.ry = 5.2426*10^6;

% ---------- Read in DEM
xyz = readBathy2vars([dem_fol fname],lim);

% ---------- Location of KKL - Load and Convert to UTM if necessary
kml_fol = 'C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\KML\RunupTransects\OwenBeach\';
kml_nm = 'OB1.kml';

temp = kml2struct([kml_fol kml_nm]); % load in KML of Transect
T.lat = temp.Lat;
T.lon = temp.Lon;
[T.x_utm,T.y_utm] = deg2utm(T.lat,T.lon); % Convert to utm
T.name = kml_nm;
clear kml_fol kmls kml_nm

% ---------- Grab x,y and depths and put in a variable
x = xyz(:,1);
y = xyz(:,2);
z = xyz(:,3);

% ---------- Name of Folder to house bathy transects
transect_folder = 'Owen_Beach_TransectBathy';

% ---------- Create bathy-transects for each transect of KML

% Set resolution of transect
ds = 1; %[m]

% Create transect
[ T.line_x, T.line_y ] = createTransect(T.x_utm(1),T.y_utm(1),T.x_utm(2),T.y_utm(2),ds);

% Extract subset of depths that covers transect
t_inds = x >= min(T.line_x) & x <= max(T.line_x) & ...
    y >= min(T.line_y) & y <= max(T.line_y);
t_x = x(t_inds);
t_y = y(t_inds);
t_z = z(t_inds);

% Find nearest depth for each location on transect
tic
T.line_z = zeros(size(T.line_x)); % DEM Line
for jj = 1:length(T.line_x)
    dist = (T.line_x(jj)-t_x).^2 + (T.line_y(jj)-t_y).^2;
    [~, I] = min(dist);
    T.line_z(jj) = t_z(I);
end
toc
clear I jj

% Reverse transect if needed, want ocean to be first
if T.line_z(1) > T.line_z(2)
    T.line_x = fliplr(T.line_x);
    T.line_y = fliplr(T.line_y);
    T.line_z = fliplr(T.line_z);
end

% Smooth the elevation data along transect
T.line_z = smoothn(T.line_z);

% Make Along Transect Grid
T.s = sqrt(T.line_x.^2+T.line_x.^2);
T.s = T.s - min(T.s);
clear t_x t_y t_z

% Save Transect
[~,temp_name,~] = fileparts(T.name);
fname = sprintf('%s_transectDepths.mat',temp_name);
temp.lat = T.lat;
temp.lon = T.lon;
temp.x_utm = T.x_utm;
temp.y_utm = T.y_utm;
temp.line_x = T.line_x;
temp.line_y = T.line_y;
temp.line_z = T.line_z;
temp.s = T.s;
temp.name = T.name;
save(fname,'-struct','temp')

% Make folder to house Transects and move transects to folder
if ~exist(transect_folder)
    mkdir(transect_folder)
    movefile(fname, transect_folder)
else
    movefile(fname, transect_folder)
end





