%RUNUP CALCULATIONS CODE
% Performs parametized Stockdon et al (2006) and TAW Runup Calculations %

% Allocate Slope variable for cat 6 we will use composite beach slope
if cat(aa) == 6
    S = S_comp(aa);
else
    S = S_beach(aa);
end

%% General parameters
% Deepwater wave length [m]
L_STK = (g*(PSP.^2))/(2*pi);      

% Mean wave period
Tm = PSP/1.1;                           % Eqn D.4.5-26

% Deep water wave lenght for the DIM method to be applied in the presence
% of structures/barriers.
L = (g*(Tm.^2))/(2*pi);    

%% Static Setup, Infragravity Swash, Incident Swash DWL and Hmo (**STK**)
% Apply the Stockdon et al. (2006) method for computing Static Setup
% and Infragravity Swash.
% Static Setup, Stockdon (2006) Eqn: 10
STK_setup = 0.35*S*sqrt(HS.*L_STK);    
 % Infragravity Swash, Stockdon (2006) Eqn: 12
STK_IG = 0.06*sqrt(HS.*L_STK);   
% Incident Swash, Stockdon (2006) Eqn: 11
STK_INC = 0.75*S*sqrt(HS.*L_STK);
% Stockdon (2006) Eqn: 7, 11, and 12
STK_swash = sqrt(STK_IG.^2 + STK_INC.^2)/2;
% Runup Stockdon (2006) Eqn: 9
STK_Runup = 1.1*(STK_setup + STK_swash);
% Total Water Levels
STK_TWL = STK_Runup + Tide;
% Determine the DWL2%
STK_DWL2 = Tide + 1.1*(STK_setup + STK_IG/2) - Ej(aa);  
% Wave heigth at toe calculated using a breaker index of 0.78
STK_Hmo = STK_DWL2*0.78;    
% If the depth limited breaking is larger than the offshore conditions,
% then the latter will be used.
STK_Hmo(STK_Hmo>HS) = HS(STK_Hmo>HS); % v0.1.4

%If the wave height at the toe is less than 0 then NaN out (NTC 16 Sep
%2013)
STK_Hmo(STK_Hmo<0) = NaN;


%% SETUP
% v0.1.5

setup = STK_setup;
DWL2 = STK_DWL2;
Hmo = STK_Hmo;

%% FLAGS
% Flag the values where STK_DWL2 is lower than Dc. For these cases we will
% use the Stockdon approach to compute the TWL. 
flag_Hmo = DWL2 <= 0;
% flag_Hmo = (Tide + 1.1*(STK_setup) - Ej(aa)) <= 0;
if sum(flag_Hmo)>0
    display(['Warning: ',num2str(sum(flag_Hmo)),...
        ' negative DWL2 values']);
    Hmo(flag_Hmo) = NaN;
end
% Flooding flag
% Check if the profile is under flooding conditions and flag it 
flag_flood = (1.1*STK_setup + Tide) >= Dc(aa);

if sum(flag_flood)>0
    display(['Warning: ',num2str(sum(flag_flood)),...
        ' conditions flood the profile (Dhigh < Tide+setup)']);
end

% If a certain category we will use a certain approach
if cat(aa) == 1
    approach = 3;
elseif cat(aa) == 2 || cat(aa) == 3
        % Use Local/TAW Combination

    approach = 1;
    
elseif cat(aa) == 6
     approach = 2;
else
    approach = 3;
    display('Warning: Unrecognized category');
    display('Warning: Stockdon et al (2006) will be used for TWLs');   
end

% If a certain category we will use a certain approach
if approach == 3
    % Use Stockdon et al (2006) for TWL
    TWL = STK_TWL;
    S_over(1:length(TWL),1) = S;
    TWL_flag(1:length(TWL),1) = 3;
elseif approach == 1 
        % Use Local/TAW Combination

    run runup_calculations_TAW_local
     TWL = TAW_TWL_local;
    S_over = S_TAW_local;
    TWL_flag(1:length(TWL),1) = 1;
    
    % Cases where the toe of the structure is not exceeded by the still water
        % level will use Stockdon et al (2006).
        TWL(flag_Hmo) = STK_TWL(flag_Hmo);
        S_over(flag_Hmo) = S;
        TWL_flag(flag_Hmo) = 3;
        S_TAW_local_flag(flag_Hmo) = NaN;
elseif approach == 2
      % Use Snsh/TAW Combination  
        run runup_calculations_TAW_snsh

        TWL = TAW_TWL_Snsh;
        S_over = S_TAW_Snsh_2;
        TWL_flag(1:length(TWL),1) = 2;
      % Unknown approach warning and go forward with Stockdon
      
      % Cases where the toe of the structure is not exceeded by the still water
        % level will use Stockdon et al (2006).
        TWL(flag_Hmo) = STK_TWL(flag_Hmo);
        S_over(flag_Hmo) = S;
        TWL_flag(flag_Hmo) = 3;
        S_TAW_local_flag(flag_Hmo) = NaN;
else
    TWL = STK_TWL;
    TWL_flag(1:length(TWL),1) = 3;
    display('Warning: Unrecognized category');
    display('Warning: Stockdon et al (2006) will be used for TWLs');   
end


% Flood cases will use Stockdon et al (2006) with beach slope depending on
% the category (see first cell block of this code) only for approach(aa)=1.
if approach==1
    TWL(flag_flood) = STK_TWL(flag_flood);
    S_over(flag_flood) = S;
    TWL_flag(flag_flood,1) = 3;
    S_TAW_local_flag(flag_flood) = NaN;
end



% Mask out the values flagged by flag_DWL2. This happens when DWL2 is lower
% than the structure toe. 
DWL2(flag_Hmo) = NaN;
Hmo(flag_Hmo) = NaN;

% For extreme water level analyses when TAW_TWL < Tide + 1.1(STK_setup +
% STK_IG/2) it will be replaced by STK_TWL
% v0.1.3

SWL_flag = TWL<(Tide+1.1*(setup+STK_IG/2));    

% for  batch processing.
DWL2actual = DWL2 + Ej(aa);                 % v0.1.3 (JA, GGM)
