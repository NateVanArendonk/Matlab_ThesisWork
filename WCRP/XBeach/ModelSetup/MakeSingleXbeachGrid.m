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


% Load Transect Data and plot - User choose elevation for dry points
T = load([transect_folder trans.name]);
% MHHW = 2.8590;
% MSL = 1.350;
% MLLW = -0.73;
% figure(1)
% plot(T.s,T.line_z)
% hold on 
% l1 = line([0 1000],[MHHW MHHW],'Color','r');
% l2 = line([0 1000],[MSL MSL],'Color','b');
% l3 = line([0 1000],[MLLW MLLW],'Color','k');
% legend([l1 l2 l3],'MHHW','MSL','MLLW')
% xlim([0 100])
% prompt = 'Please choose a dry elevation from the profile: ';
% dry_elevation = input(prompt);
% close(gcf)

% SLR Scenario?
slr = 1;
slr_val = 2.8; %{ft]

dry_elevation = 3.73;
if slr 
    slr_val = slr_val * 0.3048;
    dry_elevation = dry_elevation + slr;
end

% Reverse Transect so we have ocean first
if T.s(1) > T.s(2)
    T.s = fliplr(T.s);
end

% Find where elevation should be dry (Z > mhhw)
ind = find(T.line_z > dry_elevation);

% --------------------- Inputs for Xbeach Grid ------------------------
% Generates our grid spacing in x
xin = T.s; % Vector with cross-shore coordinates; increasing from zero - Likely the Transect
zin = T.z_s; % Vector with bed levels; positive up
tm_val = 2.34; % Incident short wave period (used to determine grid size at offshore boundary [s]
nonh_val = 1; % Creates a non-hydrostatic model grid
ppwl_val = 50; % Number of points per wave length - around 30 if domain is 1-2 wave lengths, otherwise use 100 for large domain
wl_val = 3.7231; % Water level elevation relative to bathymetry used to estimate depth [m]
if slr
    wl_val = wl_val + slr_val;
end
dx_inc = .15; % Minimum required cross shore grid size (usually over land) [m]
dry_size = 2; % Grid size for dry points [m]
dry_dist = ind(1); % Cross-shore distance from which points are dry [m]

% Inds where to subset grid for runway
depth_val = wl_val - 4; % 4 is just an aribitrary number based off of testing of toy model
inds = zin <= depth_val;

To = tm_val./1.1; % Equation 2 from Vandermeer and guidance from Robert 

% Make a runway where you make all the depth the same up until a certain
% point in the domain
make_runway = 1;
if make_runway
    zin(inds) = depth_val;
else % subset grid to be smaller domain if we aren't making a runway
    zin = zin(~inds);
    xin = xin(~inds);
    step = diff(xin); step = step(1);
    xin = 0:step:step*(length(xin)-1);  
end

% Decrease the size of the runway
limit_runway = 1;
if limit_runway 
    
    % Find Wave Length of incoming waves 
    h = 4; %meters, deep water (offshore)
    fo = 1./tm_val; % frequency
    k = getk(fo,h);
    Lo = 2*pi./k;
    
    % Lets do 2.5 wave lengths worth of runway 
    runwayLength = round(Lo*2.5);
    runway = find(inds == 1);
    if length(runway) >= runwayLength % if the runway is longer than the deisred length we can sub sample 
        r_start = runway(end) - runwayLength; % Find what indice to cut it off at 
        inds = r_start:length(xin); % Grab the indices we want
        xin = xin(inds); %sub sample
        zin = zin(inds);
        xin = 1:length(xin); % Go in and start it from 1
    end
%     inds = find(xin >= 30);
%     step = diff(xin); step = step(1);
%     xin = 0:step:step*(length(xin)-1); % not sure if -4 will be needed for all cases
%     zin = zin(inds);
end
    
% Generate Xbeach X grid
[T.sgr, T.zgr] = xb_grid_xgrid(xin, zin, 'Tm', tm_val,...
    'nonh', nonh_val, 'ppwl', ppwl_val, 'wl', wl_val,'dxmin', dx_inc,...
    'dxdry', dry_size, 'xdry',dry_dist);

% Convert s to x,y
xgr = interp1(T.s,T.line_x,T.sgr);
ygr = interp1(T.s,T.line_y,T.sgr);
zgr = T.zgr;

% Save it up 
fname = sprintf('%s_Xbeach',T.Name);
save([fname '_x.grd'],'xgr','-ascii')
save([fname '_y.grd'],'ygr','-ascii')
save([fname '_z.grd'],'zgr','-ascii')
fname_x = [fname '_x.grd'];
fname_y = [fname '_y.grd'];
fname_z = [fname '_z.grd'];

% Write delft3d xy-coordinates of grid
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


% Write JONSWAP file
% Hsig = 0.5335;
% jpath = 'C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\Tacoma\Owen\XBeach\OB_Runs\run1\';
% writeJONSWAP(Hsig,tm_val,jpath)

% input for params file  
gname = 'OB1_Xbeach_';
grdX = [gname 'x.grd'];
grdY = [gname 'y.grd'];
grdZ = [gname 'z.grd'];
grdXY = 'xy.grd';
grdLoc = 'C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\Tacoma\Owen\XBeach\Owen_Beach_XbeachGrids\';
tstop = 1200; % Seconds to run model - 3600 = 1 hour 

% Write Params file
Hsig = 0.1524;
num_runs = 'single';
parameter_testing = 'roughness';
bed_coef = 'manning';
switch num_runs
    case 'single'  % -------------- Create Single Xbeach Runs -------------
        switch bed_coef
            case 'chezy' % ----------- Chezy Roughness --------------------
                chezy_rough = 25;
                run_name = sprintf('chezy_run');
                run_fol = ['C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\Tacoma\Owen\XBeach\OB_Runs\chezy\' run_name];
                writeParams1D(grdX,grdY,grdZ,grdXY,grdLoc,tstop,wl_val,run_fol,'chezy',chezy_rough);
                % Write Jonswap
                writeJONSWAP(Hsig,tm_val,run_fol,45.0);
                
            case 'manning' % ---------- Manning Roughness -----------------
                manning_rough = .025;
                if slr
                    run_name = sprintf('man%d_SLR%.2f',manning_rough*1000,slr_val);
                else
                    run_name = sprintf('man%d_contemp',manning_rough*1000);
                end
                run_fol = ['C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\Tacoma\Owen\XBeach\OB_Runs\manning\' run_name];
                writeParams1D(grdX,grdY,grdZ,grdXY,grdLoc,tstop,wl_val,run_fol,'manning',manning_rough);
                % Write Jonswap
                writeJONSWAP(Hsig,tm_val,run_fol,45.0);
        end
        % Copy over executeables
        copyfile('C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\Tacoma\Owen\XBeach\OB_Runs\Executeables',run_fol)
        % Run Xbeach
        run = 1;
        if run
            run_xb = 'call "xbeach.exe"';
            
            cd(run_fol)
            system(run_xb)
        end
        
    case 'multiple' % ------------- Create Multiple Xbeach Runs -----------
        switch parameter_testing
            case 'roughness' % ------- Testing Influence of Roughness -----
                switch bed_coef
                    case 'chezy' % ------- Chezy Roughness ----------------
                        chezy_rough = 10:10:90;
                        for ii = 1:length(chezy_rough)
                            chezy_run = sprintf('r%d',chezy_rough(ii));
                            run_fol = ['C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\Tacoma\Owen\XBeach\OB_Runs\RoughnessTesting\chezy\' chezy_run];
                            writeParams1D(grdX,grdY,grdZ,grdXY,grdLoc,tstop,wl_val,run_fol,'chezy',chezy_rough(ii));
                            % Write Jonswap
                            writeJONSWAP(Hsig,tm_val,run_fol,45.0);
                        end
                        
                    case 'manning' % ---- Manning Roughness ---------------
                        manning_rough = 0.01:0.01:0.09;
                        for ii = 1:length(manning_rough)
                            run_name = sprintf('r%d',manning_rough(ii)*100);
                            run_fol = ['C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\Tacoma\Owen\XBeach\OB_Runs\RoughnessTesting\manning\' run_name];
                            writeParams1D(grdX,grdY,grdZ,grdXY,grdLoc,tstop,wl_val,run_fol,'manning',manning_rough(ii));
                            % Write Jonswap
                            writeJONSWAP(Hsig,tm_val,run_fol,45.0);
                        end
                end
        end
end


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
