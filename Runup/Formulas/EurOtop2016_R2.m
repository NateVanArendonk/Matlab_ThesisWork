function R2 = EurOtop2016_R2(Hs,To,beta,dlow,swl)
% Calculates run up from EuroTop II, 2002: "Manual on wave overtopping of sea defences and related structures" - equation 5a
% Hs is Hsig in deep water
% To is peak period which will be converted to zeroth moment period 
% Beta is the slope 
% dlow is the elevation of the toe of the structure/run up zone
% swl is the still water level, tid elevel at that time step
% Note Table 5.1 in manual gives guidance for choosing the correct formula - % Mean value approach vs design approach  

% Constants
g = 9.81; % m/s2
yf = 0.7; % roughness
yb = 1; % influence of berm 

% convert Peak period to zeroth moment period 
To = To/1.1; % Equation 2

% Calculate Hmo based on guidance in Tillamook study 
h = 5000; %meters, deep water (offshore)
fo = 1/To; % frequency
k = getk(fo,h); 
Lo = 2*pi/k;
n1 = 0.35*beta*sqrt(Hs*Lo);
n2 = 0.06*sqrt(Hs*Lo);
dwl = swl + 1.1*(n1 + (n2/2) - dlow);
Hmo = dwl * .78;

% Calculate Iribarren number 
so = 2*pi*Hmo/(g*To^2); % Wave steepness
Eo = beta/sqrt(so);

inds = Eo <= 1.8; % 1.8 chosen from VanderMeer figure 3
z2 = zeros(size(Hs));
z2(inds) = Hmo*(1.75*yf*Eo);
zs(~inds) = Hmo*(1.07*yf*(4.3-1.5/sqrt(yb*Eo)));
R2 = z2;
end

