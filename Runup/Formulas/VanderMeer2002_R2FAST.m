function R2 = VanderMeer2002_R2FAST(Hs,To,beta,dlow,waveAngle,transectOrientation,shorelineOrientation,swl,friction,s,z)
% Calculates run up from VanderMeer, 2002: "TAW Technical Report Wave Run-up and Wave Overtopping at Dikes" - equation 5a
% Hs is Hsig in deep water
% To is peak period which will be converted to zeroth moment period 
% Beta is the slope 
% dlow is the elevation of the toe of the structure/run up zone
% swl is the still water level, tide level at that time step 
% s is the along transect coordiantes - high resolution 
% z is the ZED level of transect - high resolution becuase we will be
% calculating slope over it if the DEM we have is good
% wave angle is angle of wave approaching shore
% shoreline orientation is general heading of shoreline, in terms of
% 'towards the' not 'from the'

% Constants
g = 9.81; % m/s2
y = friction;

% convert Peak period to zeroth moment period 
To = To./1.1; % Equation 2

% Calculate Hmo based on guidance in Tillamook study 
h = 20; %meters, deep water (offshore)
fo = 1./To; % frequency
k = getk(fo,h); 
Lo = 2*pi./k;

% !!!!!! Calculate avg slope between +- 1.5 Hs of swl - Should already be done
% prior to calling SCRIPT!!!!!
% slopeInds = find(z >= swl - Hs & z <= swl + Hs);
% if ~isempty(slopeInds)
%     % Guidance section 2.3, avg slope
%     beta = abs((z(slopeInds(end))-z(slopeInds(1)))/(s(slopeInds(end)) - s(slopeInds(1))));
% else
%     beta = .25; % Guess - relatively representative slope along ruston 
% end


n1 = 0.35.*beta.*sqrt(Hs.*Lo); % Setup from stockdon
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
HmoInds = Hmo == 0 | Hmo <= 1*10^-4;
Hmo(HmoInds) = 0;

% Calculate wave angle of attack 
shore = shorelineOrientation;
tran = transectOrientation;

% First calculate difference between transect orientation and wave angle 
window1 = [wrapTo360(tran-80);wrapTo360(tran+80)];
window2 = [wrapTo360(tran-110);wrapTo360(tran+110)];
inWindow1 = waveAngle >= window1(1) | waveAngle <= window1(2); % all waves coming within 80 degrees
inWindow2 = waveAngle >= window2(1) | waveAngle <= window2(2);
inWindow2(inWindow1) = 0; inWindow2 = logical(inWindow2); % All waves coming between 80 and 110 degrees
inWindow3 = inWindow2;
inWindow3(inWindow1) = 1;
inWindow3 = logical(inWindow3);
inWindow3 = ~inWindow3; % All waves coming greater than 110 degress 
waveAttack = NaN(length(Hs),1); % populate with wave angles 
% Find difference between transect and wave angle for first group between 0
% and 80
group1 = waveAngle;
group1(~inWindow1) = NaN; % redundant but i'll keep it 
smaller = group1<tran | group1 > tran + 130;
larger = group1>=tran & group1 < tran + 130;
group1(smaller) = wrapTo360(tran - group1(smaller));
group1(larger) = abs(tran - group1(larger));
% Find difference between transcet and wave angle for second group (80 and
% 110)
group2 = waveAngle;
group2(~inWindow2) = NaN; % redundant 
smaller = group2<tran | group2 > tran + 130;
larger = group2 >= tran & group2 < tran + 130;
group2(smaller) = wrapTo360(tran - group2(smaller));
group2(larger) = abs(tran - group2(larger));

waveAttack(inWindow1) = group1(inWindow1);
waveAttack(inWindow2) = group2(inWindow2);
% Last Group
group3 = isnan(waveAttack); % remaining group that will have wave heights be 0.1
waveAttack(group3) = 1; % Meaning no reduction factor because I'm just dropping the waves down to 10 cm
Hmo(group3) = 0.1;

% Find recution factor (Vandermeer Equation 8)
waveReduction = waveAttack;
waveReduction(inWindow1) = 1 - 0.0022*waveAttack(inWindow1);
waveReduction(inWindow2) = 1 - 0.0022*0.8;



% Calculate Iribarren number 
so = 2*pi*Hmo./(g*To.^2); % Wave steepness

% Calculate Iribarren number
Eo = beta./sqrt(so);

% Calculate R2
inds = Eo <= 1.8; % 1.8 chosen from VanderMeer figure 3
z2 = zeros(size(Hs));
z2(inds) = Hmo(inds).*waveReduction(inds).*(1.75*y*Eo(inds));
z2(~inds) = Hmo(~inds).*waveReduction(~inds).*(y*(4.3-1.6./sqrt(Eo(~inds))));
R2 = z2;
end

