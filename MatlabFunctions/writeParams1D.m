function writeParams1D(gridx,gridy,gridz,xygrid,grid_loc,tstop,wl_start,run_fol,fricCoeff,fricValue)
%Writes Params.txt file necessary for Xbeach
%   gridx = string of name of xbeach x_grid file
%   gridy = string of name of xbeach y_grid file
%   gridz = string of name of xbeach z_grid file 
%   xygrid = string of name of .grd file for xbeach 
%   grid_loc = string of path to all the above files 
%   tstop = time to stop model given as seconds
%   wl_start = water level to be prescribed on border for model run [m]
%   fricCoeff = Bed Roughness parameter - manning or chezy 
%   fricValue = Friction coefficient
fricCoeff = lower(fricCoeff);
fname = 'params.txt';
fid = fopen(fname,'w');
d = datestr(datetime('now'));

% ----------------------------- HEADER ------------------------------------
fprintf(fid,'################################################################################\n');
fprintf(fid,'### XBeach parameter settings input file                                     ###\n');
fprintf(fid,'###                                                                          ###\n');
fprintf(fid,'### date:     %s                                           ###\n',d);
fprintf(fid,'################################################################################\n');

% ---------------------------- Spacer -------------------------------------
fprintf(fid,'\n');
fprintf(fid,'\n');

% --------------------------- Model Physics -------------------------------
fprintf(fid,'### Model Physics ##############################################################\n');
fprintf(fid,'swave        = 0\n');
fprintf(fid,'flow         = 1\n');
fprintf(fid,'nonh         = 1\n');
fprintf(fid,'sedtrans     = 0\n');
fprintf(fid,'morphology   = 0\n');
fprintf(fid,'nonhq3d      = 1\n');
fprintf(fid,'################################################################################\n');

% ---------------------------- Spacer -------------------------------------
fprintf(fid,'\n');
fprintf(fid,'\n');

% -------------------------- Bed Composition ------------------------------
fprintf(fid,'### Bed composition parameters #################################################\n');
fprintf(fid,'D50          = 0.000200\n');
fprintf(fid,'D90          = 0.000300\n');
fprintf(fid,'################################################################################\n');

% ---------------------------- Spacer -------------------------------------
fprintf(fid,'\n');
fprintf(fid,'\n');

% -------------------- Flow Boundary Conditions ---------------------------
fprintf(fid,'### Flow boundary condition parameters  ########################################\n');
fprintf(fid,'front        = nonh_1d\n'); % Boundary conditions for nonhydrostatic 
fprintf(fid,'back         = abs_1d\n'); % Allows for time-varying WL to be specified at the boundary while allowing any waves propagating perpendicularly towards the boundary to be absorbed - XBEACH Manual
fprintf(fid,'################################################################################\n');

% ---------------------------- Spacer -------------------------------------
fprintf(fid,'\n');
fprintf(fid,'\n');

% ------------------------- Flow Parameters -------------------------------
switch fricCoeff
    case 'chezy' % go from 50 - 70
        fprintf(fid,'### Flow parameters  ###########################################################\n');
        fprintf(fid,'bedfriction  = chezy\n'); % Bed Friction
        fprintf(fid,'bedfriccoef  = %d\n',fricValue); 
        fprintf(fid,'################################################################################\n');
    case 'manning' % 0.02 is sandy, 0.04 is really rough 
        fprintf(fid,'### Flow parameters  ###########################################################\n');
        fprintf(fid,'bedfriction  = manning\n'); % Bed Friction
        fprintf(fid,'bedfriccoef  = %.3f\n',fricValue); 
        fprintf(fid,'################################################################################\n');
end

% ---------------------------- Spacer -------------------------------------
fprintf(fid,'\n');
fprintf(fid,'\n');


% ------------------------- Grid Parameters -------------------------------
fprintf(fid,'### Grid parameters  ###########################################################\n');
% Get Names of Grid files 
% depfile = strcat(grid_loc,'/',gridz);
% xfile = strcat(grid_loc,'/',gridx);
% yfile = strcat(grid_loc,'/',gridy);

depfile = gridz;
xfile = gridx;
yfile = gridy;
xyfile = strcat(grid_loc,'/',xygrid);

% Get Size of Grid - Note want it to be one less than grid size 
[nx,~] = getXYgridParams(xyfile);

% Write text to file 
fprintf(fid,'depfile      = %s\n',depfile);
fprintf(fid,'posdwn       = 0\n');
fprintf(fid,'alfa         0\n');
fprintf(fid,'thetamin     200\n');
fprintf(fid,'thetamax     300\n');
fprintf(fid,'dtheta       20\n');
fprintf(fid,'thetanaut    1\n');
fprintf(fid,'gridform     = xbeach\n');
fprintf(fid,'xfile        = %s\n',xfile);
fprintf(fid,'yfile        = %s\n',yfile);
fprintf(fid,'vardx        = 1\n');
fprintf(fid,'nx           = %d\n',nx);
fprintf(fid,'ny           = 0\n');
fprintf(fid,'################################################################################\n');

% ---------------------------- Spacer -------------------------------------
fprintf(fid,'\n');
fprintf(fid,'\n');


% --------------------------- Model Time ----------------------------------
fprintf(fid,'### Model time  ################################################################\n');
fprintf(fid,'tstop        = %d\n', tstop); % Seconds 
fprintf(fid,'################################################################################\n');

% ---------------------------- Spacer -------------------------------------
fprintf(fid,'\n');
fprintf(fid,'\n');


% ------------------------- Morphology Parameters -------------------------
fprintf(fid,'### Morphology parameters  #####################################################\n');
fprintf(fid,'morfac       1\n');
fprintf(fid,'wetslp       0.150000\n');
fprintf(fid,'dryslp       1\n');
fprintf(fid,'struct       0\n');
fprintf(fid,'ne_layer     nebed.dep\n');
fprintf(fid,'################################################################################\n');

% ---------------------------- Spacer -------------------------------------
fprintf(fid,'\n');
fprintf(fid,'\n');


% -------------------- Tide Boundary Conditions ---------------------------
fprintf(fid,'### Tide boundary conditions  ##################################################\n');
% Note that it appears within the tide.txt file that the three values it
% wants are time, wl and wl so we will prescribe a start and end time, and
% then the water level at each time and each point (offshore or onshore
% boundary).

% Write tide.txt file 
writeTideXbeach(tstop,wl_start,run_fol)
% Write text to params file
fprintf(fid,'zs0file      = tide.txt\n');
fprintf(fid,'tideloc      = 2\n'); % Two time varying water level signals
fprintf(fid,'################################################################################\n');

% ---------------------------- Spacer -------------------------------------
fprintf(fid,'\n');
fprintf(fid,'\n');


% -------------------- Wave Boundary Conditions ---------------------------
fprintf(fid,'### Wave boundary condition parameters  ########################################\n');
fprintf(fid,'instat       = jons\n');
fprintf(fid,'wbcversion   = 3\n');
fprintf(fid,'bcfile       = jonswap\n');
fprintf(fid,'random       = 0\n');
fprintf(fid,'rt           = %d\n',tstop); % Duration of wave spectrum at off-shore boundary, in morphological time [s]
fprintf(fid,'dtbc         1\n');
% Random turns on/off random seed
% rt can be less than tsop, but for run-up statistics, no better
% Trick is to use time series for jonswap spec and turn random on,
% Hilbert transform is exponential in computation, rt is 600 is 1/4 time rt is 1200
fprintf(fid,'################################################################################\n');

% ---------------------------- Spacer -------------------------------------
fprintf(fid,'\n');
fprintf(fid,'\n');


%%% Output variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tstart should be greater than zero, let the model "burn in"
% for run up stats want to sample at smallish fraction of wave period
%zb		% Bed level
%zs		% Water Level
%u		% Glm velocity in cell centre, x-component
%qx		% Discharge in u-points, x-component
%H 		% Hrms wave height based on instantaneous wave energy
%Hrunup 	% Short wave height used in runup formulation
%Qb 		% Fraction breaking waves
%hh 		% Water Depth
%taubx 		% Bed shear stress, x-component
%tauby 		% Bed shear stress, y-component
%D 		% Dissipation
%E 		% Wave Energy
%L1 		% Wave Length used in dispersion relation
%iwl 		% Location of water line (including long wave runup
%kb 		% Near bed turbulence intensity due to depth induces breaking
%runup 		% Short wave runup height
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(fid,'### Output  ####################################################################\n');
fprintf(fid,'outputformat = netcdf\n');
fprintf(fid,'tintm        = 300\n');
fprintf(fid,'tintg        = 1\n');
fprintf(fid,'tstart       = 0\n');
fprintf(fid,'tintp	      = 0.25\n');
fprintf(fid,'\n');
fprintf(fid,'nglobalvar = 2\n');
fprintf(fid,'zb\n');
fprintf(fid,'zs\n');
fprintf(fid,'\n');
fprintf(fid,'nrugauge = 1\n');
fprintf(fid,'0 0\n'); % Start gauge anywhere
fprintf(fid,'################################################################################\n');
% Close Params File
fclose(fid);
fclose('all');

% Move file 
movefile('params.txt',run_fol)
end

