clearvars
addpath C:\Users\ahooshmand\Desktop\Xbeach\XBeach_Codes
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
kmls = dir('C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\KML\RunupTransects\OwenBeach\*.kml');
for ii = 1:length(kmls)
    kml_nm = kmls(ii).name;
    temp = kml2struct([kml_fol kml_nm]); % load in KML of Transect
    T(ii).lat = temp.Lat;
    T(ii).lon = temp.Lon;
    [T(ii).x_utm,T(ii).y_utm] = deg2utm(T(ii).lat,T(ii).lon); % Convert to utm
    T(ii).name = kml_nm;
end
clear kml_fol kmls kml_nm

% ---------- Grab x,y and depths and put in a variable 
x = xyz(:,1);
y = xyz(:,2);
z = xyz(:,3);

% ---------- Name of Folder to house bathy transects
transect_folder = 'Owen_Beach_TransectBathy';

% ---------- Create bathy-transects for each transect of KML
for ii = 1:length(T)
    
    % Set resolution of transect
    ds = 1; %[m]
    
    % Create transect
    [ T(ii).line_x, T(ii).line_y ] = createTransect(T(ii).x_utm(1),T(ii).y_utm(1),T(ii).x_utm(2),T(ii).y_utm(2),ds);

    % Extract subset of depths that covers transect
    t_inds = x >= min(T(ii).line_x) & x <= max(T(ii).line_x) & ...
        y >= min(T(ii).line_y) & y <= max(T(ii).line_y);
    t_x = x(t_inds);
    t_y = y(t_inds);
    t_z = z(t_inds);
 
    % Find nearest depth for each location on transect
    tic
    T(ii).line_z = zeros(size(T(ii).line_x)); % DEM Line
    for jj = 1:length(T(ii).line_x)
        dist = (T(ii).line_x(jj)-t_x).^2 + (T(ii).line_y(jj)-t_y).^2;
        [~, I] = min(dist);
        T(ii).line_z(jj) = t_z(I);
    end
    toc
    clear I jj 
    
    % Reverse transect if needed, want ocean to be first
    if T(ii).line_z(1) > T(ii).line_z(2)
        T(ii).line_x = fliplr(T(ii).line_x);
        T(ii).line_y = fliplr(T(ii).line_y);
        T(ii).line_z = fliplr(T(ii).line_z);
    end
    
    % Smooth the elevation data along transect
    T(ii).line_z = smoothn(T(ii).line_z);
    
    % Make Along Transect Grid
    T(ii).s = sqrt(T(ii).line_x.^2+T(ii).line_x.^2);
    T(ii).s = T(ii).s - min(T(ii).s);    
    clear t_x t_y t_z
    
    % Save Transect
    [~,temp_name,~] = fileparts(T(ii).name);
    fname = sprintf('%s_transectDepths.mat',temp_name);
    temp.lat = T(ii).lat;
    temp.lon = T(ii).lon;
    temp.x_utm = T(ii).x_utm;
    temp.y_utm = T(ii).y_utm;
    temp.line_x = T(ii).line_x;
    temp.line_y = T(ii).line_y;
    temp.line_z = T(ii).line_z;
    temp.s = T(ii).s;
    temp.name = T(ii).name;
    save(fname,'-struct','temp')
    
    % Make folder to house Transects and move transects to folder 
    if ~exist(transect_folder)
        mkdir(transect_folder)
        movefile(fname, transect_folder)
    else
        movefile(fname, transect_folder)
    end
end





