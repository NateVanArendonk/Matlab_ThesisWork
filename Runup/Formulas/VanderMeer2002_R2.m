function R2 = VanderMeer2002_R2(Hs,To,beta,depth,dlow,swl,friction,s,z)
% Calculates run up from VanderMeer, 2002: "TAW Technical Report Wave Run-up and Wave Overtopping at Dikes" - equation 5a
% Hs is Hsig in deep water
% To is peak period which will be converted to zeroth moment period 
% Beta is the slope 
% dlow is the elevation of the toe of the structure/run up zone
% swl is the still water level, tid elevel at that time step 
% s is the along transect coordiantes - high resolution 
% z is the ZED level of transect - high resolution becuase we will be
% calculating slope over it if the DEM we have is good

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

% Calculate avg slope between +- 1.5 Hmo of swl 
slopeInds = find(z >= swl - Hmo & z <= swl + Hmo);
if ~isempty(slopeInds)
    % Guidance section 2.3, avg slope
    beta = abs((z(slopeInds(end))-z(slopeInds(1)))/(s(slopeInds(end)) - s(slopeInds(1))));
else
    beta = .3; % Guess 
end

% Calculate Iribarren number
Eo = beta./sqrt(so);

% Calculate R2
inds = Eo <= 1.8; % 1.8 chosen from VanderMeer figure 3
z2 = zeros(size(Hs));
z2(inds) = Hmo(inds).*(1.65*y*Eo(inds));
z2(~inds) = Hmo(~inds).*(y*(4.3-1.6./sqrt(Eo(~inds))));

% --------------------- Second Iteration with New Slope -------------------
% Now that R2 is known, rerun for new range with -1.5Hmo and +R2 above SWL
slopeInds = find(z >= swl - Hmo & z <= swl + z2);
if ~isempty(slopeInds)
    % Guidance section 2.3, avg slope
    beta = abs((z(slopeInds(end))-z(slopeInds(1)))/(s(slopeInds(end)) - s(slopeInds(1))));
else
    beta = .3; % Guess 
end


% Calculate Iribarren number
Eo = beta./sqrt(so);

inds = Eo <= 1.8; % 1.8 chosen from VanderMeer figure 3
z2 = zeros(size(Hs));
z2(inds) = Hmo(inds).*(1.75*y*Eo(inds));
z2(~inds) = Hmo(~inds).*(y*(4.3-1.6./sqrt(Eo(~inds))));

R2 = z2;
end

