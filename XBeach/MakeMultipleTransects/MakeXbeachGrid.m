clearvars
addpath C:\Users\ahooshmand\Desktop\Xbeach\XBeach_Codes
addpath C:\Functions_Matlab
addpath C:\Functions_Matlab\mapping\kml

% ADD OET!!!!!!
dry_elevation = 2.73; %[m]

% Get Contents of Transect folder
transect_folder = 'Owen_Beach_TransectBathy\';
trans = dir('Owen_Beach_TransectBathy\*.mat');

% Name of Grid folder 
grid_folder = 'Owen_Beach_XbeachGrids';

% Loop through transects and make grids 
for ii = 1:length(trans)
    
    % Load Transect Data
    T(ii) = load([transect_folder trans(ii).name]);
    
    % Reverse Transect so we have ocean first 
    if T(ii).s(1) > T(ii).s(2)
        T(ii).s = fliplr(T(ii).s);
    end
    
    % Find where elevation should be dry (Z > mhhw)
    ind = find(T(ii).line_z > dry_elevation); 
    
    % --------------------- Inputs for Xbeach Grid ------------------------
    % Generates our grid spacing in x
    xin = T(ii).s; % Vector with cross-shore coordinates; increasing from zero - Likely the Transect
    zin = T(ii).line_z; % Vector with bed levels; positive up
    tm_val = 4/1.1; % Incident short wave period (used to determine grid size at offshore boundary [s]
    nonh_val = 1; % Creates a non-hydrostatic model grid
    ppwl_val = 50; % Number of points per wave length - around 30 if domain is 1-2 wave lengths, otherwise use 100 for large domain
    wl_val = 1; % Water level elevation relative to bathymetry used to estimate depth [m]
    dx_inc = .15; % Minimum required cross shore grid size (usually over land) [m]
    dry_size = 2; % Grid size for dry points [m]
    dry_dist = ind(1); % Cross-shore distance from which points are dry [m]
    
    % Generate Xbeach X grid
    [T(ii).sgr, T(ii).zgr] = xb_grid_xgrid(xin, zin, 'Tm', tm_val,...
'nonh', nonh_val, 'ppwl', ppwl_val, 'wl', wl_val,'dxmin', dx_inc,...
'dxdry', dry_size, 'xdry',dry_dist);
    
    % Convert s to x,y
    xgr = interp1(T(ii).s,T(ii).line_x,T(ii).sgr);
    ygr = interp1(T(ii).s,T(ii).line_y,T(ii).sgr);
    zgr = T(ii).zgr;
    
    fname = sprintf('%s_Xbeach',T(ii).Name);
    save([fname '_x.grd'],'xgr','-ascii')
    save([fname '_y.grd'],'ygr','-ascii')
    save([fname '_z.grd'],'zgr','-ascii')
    fname_x = [fname '_x.grd'];
    fname_y = [fname '_y.grd'];
    fname_z = [fname '_z.grd'];
    
    % Make folder to house Transect Grids and move to folder
    if ~exist(grid_folder)
        mkdir(grid_folder)
        movefile(fname_x, grid_folder)
        movefile(fname_y, grid_folder)
        movefile(fname_z, grid_folder)
    else
        movefile(fname_x, grid_folder)
        movefile(fname_y, grid_folder)
        movefile(fname_z, grid_folder)
    end
    clear fname_x fname_y fname_z 
    
    plotting = 1;
    if plotting
        clf
        hold on
        plot(T(ii).s,T(ii).line_z) % Plot the transect
        plot(T(ii).sgr,zgr,'*') % Plot the grid points on the transect
        plot(T(ii).sgr(2:end),diff(T(ii).sgr),'ro')
        
        pause
        clf
        scatter(xgr,ygr,[],zgr)
        colorbar
    end
end
