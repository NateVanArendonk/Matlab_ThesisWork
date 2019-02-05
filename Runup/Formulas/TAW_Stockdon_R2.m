function R2 = TAW_Stockdon_R2(Ho,Tp,beta,Dlow,SWL,roughness)
%Calculation of Run-up (R2) using the TAW method described in the tillamook
% coastal hazards study

% First calculate DWL - dynamic water level
% SWL is the water height at location of run-up calculation - tide level
% Dlow is toe depth of location
% Basically the vandermeer formula with stockdon to calculate Hmo

yf = roughness;


% Calculate Lo
h = 5000; %meters, deep water (offshore)
fo = 1/Tp;
k = getk(fo,h); 
Lo = 2*pi/k;

% First calculate DWL - dynamic water level and Hmo
n_s = .35*beta*sqrt(Ho*Lo); % Setup from stockdon
n_r = 0.06*sqrt(Ho*Lo); % IG Swash from stockdon 
DWL = SWL + 1.1 * (n_s + (n_r/2)) - Dlow; 
Hmo = DWL * 0.78;

% Shoaling Factor - Ecg = Ecg
[~,Cgshallow] = getcp(fo,abs(SWL-Dlow));
[~,Cgdeep] = getcp(fo,5000);
Htoe = ((Cgdeep/Cgshallow)^2)*Ho;
% Take smaller of two to be the wave height at the toe of the structure 
Hmo = min(Htoe,Hmo);

IB = beta/sqrt(Hmo/Lo);

inds = IB <= 1.8; % 1.8 chosen from VanderMeer figure 3
R2 = zeros(size(Hmo));

% Calculate Run-up
R2(inds) = Hmo*(1.75*yf*IB);
R2(~inds) = Hmo * (1*yf*(4.3 - (1.6/sqrt(IB))));
end

