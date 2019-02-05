function R2 = VanderMeer2002_WCRPAdapted_R2(Hs,To,beta,depth,dlow,swl,friction)
% Calculates run up from VanderMeer, 2002: "TAW Technical Report Wave Run-up and Wave Overtopping at Dikes" - equation 5a
% Hs is Hsig in deep water
% To is peak period which will be converted to zeroth moment period 
% Beta is the slope 
% dlow is the elevation of the toe of the structure/run up zone
% swl is the still water level, tid elevel at that time step 

% Constants
g = 9.81; % m/s2
y = friction;

% convert Peak period to zeroth moment period 
To = To./1.1; % Equation 2

% Calculate Hmo based on guidance in Tillamook study 
h = depth; %meters, deep water (offshore)
fo = 1./To; % frequency
k = getk(fo,h); 
Lo = 2*pi./k;

n1 = 0.35*beta*sqrt(Hs.*Lo); % Setup from stockdon
n2 = 0.06*sqrt(Hs.*Lo); % IG Swash from stockdon
dwl = swl + 1.1*(n1 + (n2/2)) - dlow;
Hmo = dwl .* 0.78;

% Shoaling Factor - Ecg = Ecg
[~,Cgshallow] = getcp(fo,abs(swl-dlow));
[~,Cgdeep] = getcp(fo,5000);
Htoe = ((Cgdeep./Cgshallow).^2).*Hs;
% Take smaller of two to be the wave height at the toe of the structure 
Hmo = min(Htoe,Hmo);

% If the wave height being prescribed is zero or super fucking small, just quit here and R2 = 0
if Hmo == 0 || Hmo <= 1*10^-4
    R2 = 0;
    return
end

% Calculate Iribarren number 
so = 2*pi*Hmo./(g*To.^2); % Wave steepness

% Calculate Iribarren number
Eo = beta./sqrt(so);

% Calculate R2
inds = Eo <= 1.8; % 1.8 chosen from VanderMeer figure 3
z2 = zeros(size(Hs));
z2(inds) = Hmo(inds).*(1.65*y*Eo(inds));
z2(~inds) = Hmo(~inds).*(y*(4.3-1.6./sqrt(Eo(~inds))));
R2 = z2;
end

