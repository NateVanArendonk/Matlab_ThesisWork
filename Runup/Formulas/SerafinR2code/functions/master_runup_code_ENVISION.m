%% MASTER CODE TO PERFORM RUNUP AND OVERTOPPING CALCULATIONS %
% Clean Matlab
close all
clc
clear all

%% Initial parameter selection and calculations         (v0.1, v0.1.1)
% Constants
g = 9.81;                       % Gravity (SI) [m/s2]
g_ft = 32.174;                  % Gravity (US) [ft/s2]
m2ft = 3.281;                   % meters to feet conversion
ft2m = 0.3048;                  % feet to meters conversion
rho = 1025;                     % Density [kg/m3];
gamma = 0.78;                   % Saturated Breaking Coefficient
vdatum = 0.0;                  % Vertical datum adjustment [m] 
                                % Relative elevation of profile datum wrt
                                % NAVD88 THIS DATA I HAVE PUT IN IS
                                % RELATIVE TO MLLW SO THIS IS 0.108. IF
                                % RELATIVE TO NAVD88 vdatum = 0;
                                
% Input Wave parameters
load(['C:\data\projects\coca\Runup_Formulations\data\tillamook\',...
    'sampleWaveData.mat']);


% Input Wave parameters
load(['Z:\root\home\shusin3\shared\',...
    'dogami_fema\tillamook\wave_data\event_forcing.mat']);

% Bathymetry - All the transects to be considered should be in the same
% folder or at least soft-linked to it. They are assumed to be named
% according to the PRIDs
% Select the parent folder for the bathymetry

bathy_folder = ['C:\data\projects\coca\Runup_Formulations\data\tillamook\'...
    'transects\bathy\'];

% Prefix to identify the bathymetry files and the output files
prefix = 'TILL';

% Profile data (to be read from an MS Excel Table)
% profile_data = ['/home/server/pi/homes/ncohn/FEMA/',...
%     'lincoln/transects/02-table/',...
%     'Lincoln_insitucombined_EVT_Results_Transect_info_111213.xlsx'];
profile_data = ['C:\data\projects\coca\Runup_Formulations\data\tillamook\'...
    'Tillamook_insitucombined_EVT_Results_Transect_info_030813.xlsx'];
profile_sheet = 'Sheet1';
profile_range = 'E4:R181';



% Output Folder ~~~~~~~~~~~~~~~~~~~~~~~~~
outputfolder = ['C:\Users\kserafin\Documents\ENVISIONtesting\'];
% Folder to save the results of this code

% If you want a map with the discretized plots add the path to the
% coastline database. If no map is desired define as an empty variable.
coastline = 'C:\data\projects\coca\Runup_Formulations\data\tillamook\WCSL_UTM_FIN.mat';


%% Import data for computations
prof = [45,76,77,78,79,80,81,93,109,110,111,176,177,178];

% Get starting time
start_time = now;

% Read the profile data
Data_profile = xlsread(profile_data, profile_sheet, profile_range);
Data_profile(Data_profile==-999)=NaN;
% Assign variables (verify the Excel sheet)
cat = Data_profile(:,1);                % Get the profile category
BBC=Data_profile(:,2);                  % barrier/beach crest location - (overtopping).        
Dc=Data_profile(:,3);                   % critical 'barrier'/dune crest heights (Dhigh)
Ej=Data_profile(:,4);                   % Beach dune juncture elevation (Dlow)
S_beach = Data_profile(:,5);            % Beach slope (tan(beta))
S_comp = Data_profile(:,6);             % Composite beach slope (tan(beta))
bBackshore=Data_profile(:,7);           % backshore elevation drop to point where backshore slope angle derived -(overtopping)
mBackshore=Data_profile(:,8);           % backshore slope angle (overtopping)
mOverflow = Data_profile(:,9);         % Slope for bore propagation (overland flow), after structure has overtopped. 
SA=Data_profile(:,10);                  % shoreline angle (Wave Direction)
yr = Data_profile(:,11);                % Roughness factor 
yb_or = Data_profile(:,12);             % Berm reduction factor (overwritte value)
approach = Data_profile(:,14);          % Approach flag (Local/TAW(1), Snsh/TAW(2) or STK(3))

yb_or(:) = 1;
% Bathymetry identifiers (PRID)
prid = xlsread(profile_data,profile_sheet,'B4:B181');

%% Calculate Nearshore Slopes - Used in parameterized DIM calculations below %%

% Preallocate some variables
twl_max = nan(length(cat),1);
q_max = nan(size(twl_max));
fc_min = nan(size(twl_max));
hv2_max = nan(size(twl_max));

% Loop over transects, variable aa will be passed to the other codes.
for aa = prof(1):prof(end)
    % Close all the generated figures
    %close all
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    % Load hydrodynamic data for each transect here             (v.0.1.1)
      HS = wave_height; %wave data
      WDir = wave_dir; % direction data
      PSP = wave_period; % period data
      Tide = waterlevel; %water level data (tide +ntr)
      event_time_fin = time; %time
      
      ind = find(isnan(WDir)>0);
      HS(ind) = [];
      WDir(ind) = [];
      PSP(ind) = [];
      Tide(ind) = [];
      event_time_fin(ind) = [];
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
    
    % Compute wave properties according to LWT                  (v0.1)
    % calculate wave length, breaker height, and breaking depth
    L = (g*(PSP.^2))/(2*pi);                % Deepwater wave length (m)
    Hb = 0.39*(g^0.2)*(PSP.*HS.^2).^0.4;    % Calculate the wave breaker height (m)
    hb = Hb./gamma;                         % Calculate breaking depth (m)
    
    % Check the existence of the bathymetry file                (v0.1.1)               
    chkdir = exist([bathy_folder,prefix,num2str(prid(aa)),...
        '_bathy.mat'],'file');
    
    % Go into the loop the bathymetry file is found
    if chkdir == 2
        
        % Display the considered file                           (v0.1.1)
        display(' ');
        display(['Evaluating ',bathy_folder,prefix,num2str(prid(aa)),...
            '_bathy.mat']);
        display(['File ',num2str(aa),' of ',num2str(length(cat))]);
        
        % Import topo data                                      (v0.1.1)
        bathy = load([bathy_folder,prefix,num2str(prid(aa)),'_bathy.mat']);
        x = bathy.x;            % Allocate variables
        y = bathy.y;            % Allcoate variables
        
        % Find non data points and clear them
        nanid = isnan(y);
        y(nanid) = [];
        x(nanid) = [];            
        lwt_setup_profile = nan(length(HS), length(y)-2000);
        % Adjust for the vertical datum                         (v0.1.1)
        y = y + vdatum;       
        
        % Perform runup calculations -------------------------
        % return
        run runup_calculations_ENVISION   

        % Save all the initial data in binary (mat) format --------------
        
        % Get the code information
        information.created = ['File created at ',datestr(now),' PT'];
        information.script = [mfilename('fullpath'),'.m'];
        
        % Create data array
        hydraulic_data.HS = HS;
        hydraulic_data.PSP = PSP;
        hydraulic_data.Tm = Tm;
        hydraulic_data.L = L;
        hydraulic_data.WDir = WDir;       
        hydraulic_data.Tide = Tide;
        hydraulic_data.event_time = event_time_fin;
       
%         hydraulic_data.beachslope = S;
%         hydraulic_data.Snsh = S_TAW_Snsh_2;
%         hydraulic_data.localslope = S_TAW_local;
       
        hydraulic_data.STK_Runup = STK_Runup;       
        hydraulic_data.STK_swash = STK_swash;
        hydraulic_data.STK_setup = STK_setup;       
        hydraulic_data.STK_TWL = STK_TWL;
        hydraulic_data.STK_DWL2 = STK_DWL2;        
        hydraulic_data.STK_Hmo = STK_Hmo;        
        hydraulic_data.STK_IG = STK_IG;
        hydraulic_data.STK_INC = STK_INC;
        hydraulic_data.approach = approach;

         
        if approach == 3
            hydraulic_data.beachslope = S;
        elseif approach == 1 
             hydraulic_data.localslope = S_TAW_local;
             hydraulic_data.TAW_TWL_local = TAW_TWL_local;
             hydraulic_data.TAW_swash_local = TAW_R_local;
             hydraulic_data.S_TAW_local_flag = S_TAW_local_flag;   


        elseif approach == 2
            hydraulic_data.Snsh = S_TAW_Snsh_2;
            hydraulic_data.TAW_TWL_Snsh = TAW_TWL_Snsh;
            hydraulic_data.TAW_swash_Snsh = TAW_R_Snsh_2;
            hydraulic_data.S_TAW_Snsh_flag = S_TAW_Snsh_flag;
  

        end
        


                       
        hydraulic_data.DWL2 = DWL2;
        hydraulic_data.DWL2actual = DWL2actual;
        hydraulic_data.Hmo = Hmo;        
        hydraulic_data.TWL = TWL;       
        hydraulic_data.TWL_flag = TWL_flag;
        hydraulic_data.SWL_flag = SWL_flag;
               
        % Create variable information structure 
        information.variables = ...
            {'HS = wave height (m)';...
            'PSP = Peak Wave Period (s)';...
            'Tm = Mean wave period (s), computed as PSP/1.1';...
            'L = Deepwater Wave Length (m)';...
            'WDir = Deepwater Incoming Wave Direction (degrees) Nautical Convention';...
            'Tide = Corresponding tide (m)';...
            'event_time = time associated to the TWL time series (matlab serial time)';
            'beachslope = slope used for computations (m/m)';...
            'Snsh = Slope from TAW_TWL to [Tide - 1.5Hmo] (m/m)';...
            'localslope = Slope from STK_TWL to [1.1*STK_setup + Tide] (m/m)';...
            'STK_Runup = 2% runup from Stockdon et al. (2006) Equation 19';...
            'STK_swash = Full Swash Oscillation Stockdon (2006), sqrt(STK_IG^2 + STK_INC^2)/2';...
            'STK_setup = Static setup (etabar) from Stockdon et al. 2006, Eqn: 10';...
            'STK_TWL = Total Water Level from Stockdon et al. (2006), STK_Runup + Tide';...
            'STK_DWL2 = STK_setup + STK_swash + Tide';...
            'STK_Hmo = Wave height at toe of feature, calculated with breaker index = 0.78';...
            'STK_IG = Infragravity swash (Stockdon et al. (2006) Eqn: 12)';...
            'STK_INC = Incident swash (Stockdon et al. (2006) Eqn: 11';...
            'TAW_TWL_Snsh = Total water level with the TAW approach (m) and Snsh';...
            'TAW_swash_Snsh = Swash oscillations computed with the TAW approach (m) and Snsh';...
            'TAW_TWL_local = Total water level with the TAW approach (m) and local structure slope';...
            'TAW_swash_local = Swash oscillations computed with the TAW approach (m) and local structure slope';...
            'LWT_DWL2 = LWT_setup + STK_IG + Tide [m]';...
            'LWT_Hmo = LWT_DWL2 - Ej [m]';...
            'LWT_setup = static setup at the dune toe computed using linear wave theory [m]';...
            'DWL2 = Final 2% runup level over the dune/structure toe, based on profile category and slope';...
            'DWL2actual = Final 2% runup level, based on profile category and slope';...
            'Hmo = Final wave height at toe of feature based on profile category and slope';...
            'TWL = Final total water level (m), selection based on profile category and slope';...
            'TWL_flag = Approach flag: Local Slope (1), Snsh (2), Stockdon (3), and LWT/Local (4)';...
            'S_TAW_local_flag = Slope flag for local slope approach: Ytop = STK_TWL (1) and Ytop = Dhigh (2)';...
            'S_TAW_Snsh_flag = Slope flag for Snsh slope approach: Ytop = Tide + 1.5*STK_Hmo (1) and Ytop = Dhigh (2)';...
            'SWL_flag = if TAW_TWL < Tide+1.1(STK_setup + STK_IG/2) then TWL = STK_TWL'};
        
%         % need to write hydraulic data so they are relevant to the specific profile of interest.
%         save([out_folder,prefix,num2str(prid(aa)),'_hydraulic_data.mat'],...
%             'hydraulic_data','information');
        
        % Find the maximum TWL and Freeboard
        twl_max(aa) = nanmax(TWL);
       
      
        % Plot all TWL with respect to the dune elevation
        figure
        hold on            
        plot([x(1) x(end)],[TWL TWL])
        plot(x,y,'k','LineWidth',2);
        title(['Site ',prefix,num2str(prid(aa))]);
        hold off
               
     else
      
        % Waring message                                        (v0.1.1)
        display(' ');   
        display(['File ',num2str(aa),' of ',num2str(length(cat))]);
        display('Warning:  Bathymetry File not found');
        display(['Warning:  ', bathy_folder,prefix,...
            num2str(prid(aa)),'_bathy.mat']);
    end
   
    
end

