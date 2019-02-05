close all
clearvars
addpath C:\Functions_Matlab
addpath C:\Functions_Matlab\mapping\kml

% ADD OET!!!!!!


% Get Contents of Transect folder
transect_folder = 'Owen_Beach_TransectBathy\';
trans = dir('Owen_Beach_TransectBathy\*.mat');

% Name of Grid folder to save in to
grid_folder = 'Owen_Beach_XbeachGrids';

% Loop through transects and make grids

% Load Transect Data and plot - User choose elevation for dry points
T = load([transect_folder trans.name]);
MHHW = 2.8590;
MSL = 1.350;
MLLW = -0.73;
figure(1)
plot(T.s,T.line_z)
hold on 
l1 = line([0 1000],[MHHW MHHW],'Color','r');
l2 = line([0 1000],[MSL MSL],'Color','b');
l3 = line([0 1000],[MLLW MLLW],'Color','k');
legend([l1 l2 l3],'MHHW','MSL','MLLW')
xlim([0 100])
prompt = 'Please choose a dry elevation from the profile: ';
dry_elevation = input(prompt);
close(gcf)
% Reverse Transect so we have ocean first
if T.s(1) > T.s(2)
    T.s = fliplr(T.s);
end

% Find where elevation should be dry (Z > mhhw)
ind = find(T.line_z > dry_elevation);

% --------------------- Inputs for Xbeach Grid ------------------------
% Generates our grid spacing in x
xin = T.s; % Vector with cross-shore coordinates; increasing from zero - Likely the Transect
zin = T.line_z; % Vector with bed levels; positive up
tm_val = 4; % Incident short wave period (used to determine grid size at offshore boundary [s]
nonh_val = 1; % Creates a non-hydrostatic model grid
ppwl_val = 30; % Number of points per wave length - around 30 if domain is 1-2 wave lengths, otherwise use 100 for large domain
wl_val = 1; % Water level elevation relative to bathymetry used to estimate depth [m]
dx_inc = .15; % Minimum required cross shore grid size (usually over land) [m]
dry_size = 2; % Grid size for dry points [m]
dry_dist = ind(1); % Cross-shore distance from which points are dry [m]

% Inds where to subset grid
depth_val = 4.5;
inds = zin <= (dry_elevation - depth_val);


% Make a runway where you make all the depth the same up until a certain
% point in the domain
make_runway = 0;
if make_runway
    zin(inds) = dry_elevation - depth_val;
else % subset grid to be smaller domain if we aren't making a runway
    zin = zin(~inds);
    xin = xin(~inds);
    step = diff(xin); step = step(1);
    xin = 0:step:step*(length(xin)-1);  
end

% Decrease the size of the runway
limit_runway = 0;
if limit_runway 
    inds = find(xin >= 30);
    step = diff(xin); step = step(1);
    xin = 0:step:step*(length(xin)-1); % not sure if -4 will be needed for all cases
    zin = zin(inds);
end
    
% Generate Xbeach X grid
[T.sgr, T.zgr] = xb_grid_xgrid(xin, zin, 'Tm', tm_val,...
    'nonh', nonh_val, 'ppwl', ppwl_val, 'wl', wl_val,'dxmin', dx_inc,...
    'dxdry', dry_size, 'xdry',dry_dist);

% Convert s to x,y
xgr = interp1(T.s,T.line_x,T.sgr);
ygr = interp1(T.s,T.line_y,T.sgr);
zgr = T.zgr;

fname = sprintf('%s_Xbeach',T.Name);
save([fname '_x.grd'],'xgr','-ascii')
save([fname '_y.grd'],'ygr','-ascii')
save([fname '_z.grd'],'zgr','-ascii')
fname_x = [fname '_x.grd'];
fname_y = [fname '_y.grd'];
fname_z = [fname '_z.grd'];

% Write delft3d xy-coordinates of calculation grid
proj_name = 'xy';
wlgrid('write','Filename',proj_name,'X',xgr,'Y',ygr','AutoEnclosure')
% wldep('write','bed.dep',zgr)

% Make folder to house Transect Grids and move to folder
if ~exist(grid_folder)
    mkdir(grid_folder)
    movefile(fname_x, grid_folder)
    movefile(fname_y, grid_folder)
    movefile(fname_z, grid_folder)
    movefile('xy.enc',grid_folder)
    movefile('xy.grd',grid_folder)
else
    movefile(fname_x, grid_folder)
    movefile(fname_y, grid_folder)
    movefile(fname_z, grid_folder)
    movefile('xy.enc',grid_folder)
    movefile('xy.grd',grid_folder)
end
clear fname_x fname_y fname_z

plotting = 0;
if plotting
    clf
    hold on
    plot(T.s,T.line_z) % Plot the transect
    plot(T.sgr,zgr,'*') % Plot the grid points on the transect
    plot(T.sgr(2:end),diff(T.sgr),'ro')
    
    pause
    clf
    scatter(xgr,ygr,[],zgr)
    colorbar
end

