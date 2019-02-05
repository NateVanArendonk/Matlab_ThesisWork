function [R2, IB, setup, R2_IG] = Stockdon2006_R2( Ho, To, beta, waveAngle, transectOrientation)
% [R2, IB] = stockdonR2( Ho, To, beta )
%  Estimate R2% run-up level using Stockdon et al 2006 formulation, Eq: 19
%  Also estimates Iribarren number, IB
%
% INPUTS
%   Ho - deep water sig wave height
%   To - deep wave wave period
%   beta - Slope where beta = rise / run (positive)
%       note: beta = tan(alpha) where alpha is angle
%       Stockdon defines beta_f as +/- 2sig around mean of WL timeseries 

% RS(ii).high1 = round(Stockdon2006_R2(W(ii).R1MLErth*0.3048,W(ii).TP1,0.3)*3.28084,2);


% Convert fo to Lo
h = 15; %meters, deep water (offshore)
fo = 1./To;
k = getk(fo,h); 
Lo = 2*pi./k;
dimsLo = size(Lo);
dimsK = size(k);
if dimsLo(1) ~= dimsK(1)
    Lo = Lo';
end


% ------------ Calculate wave angle of attack -----------------------------
% Added based on guidance from TAW 
tran = transectOrientation;

% Sorry this part is super confusing and I didn't comment it enough 
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
waveAttack = NaN(length(Ho),1); % populate with wave angles 
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
group3 = isnan(waveAttack); % remaining group that will have wave heights be 0
waveAttack(group3) = 0; % Meaning no reduction factor because I'm just dropping the waves down to 10 cm

% Find recution factor (Vandermeer Equation 8)
waveReduction = waveAttack;
waveReduction(inWindow1) = 1 - 0.0022*waveAttack(inWindow1);
waveReduction(inWindow2) = 1 - 0.0022*0.8;
% -------------------------------------------------------------------------


% Estimate Iribarren number with offshore conditions
% Note Beta is the average slope over 2 sigma of the water level timeseries
IB = beta./sqrt(waveReduction.*Ho./Lo);

% Calc setup
setup = 0.35*beta.*sqrt(waveReduction.*Ho.*Lo);

% Calc IG dominated run-up where IB < 0.3
R2_IG = 0.043.*sqrt(waveReduction.*Ho.*Lo);

% Calc R2%, combines Setup and Swash
R2 = 1.1*(0.35*beta.*sqrt(waveReduction.*Ho.*Lo) + sqrt(waveReduction.*Ho.*Lo.*(0.56*beta^2 + 0.004))/2);

end


function [cp, cg] = getcp(f,h)
% [cp, cg] = getc(f,h)
g = 9.81;
k = getk(f,h);
cp = sqrt(g./k.*tanh(k.*h));
%cg = 0.5*(g*tanh(k*h)+g*(k*h)*(sech(k*h)^2))/sqrt(g*k*tanh(k*h));
% cg = cp./2*(1+2*k.*h/sinh(2*k.*h)); %Kundu & Cohen p240
cg = cp./(2*(1+2*k.*h./sinh(2*k.*h)));

end

function k = getk(f,h) 
%k = getk(f,h)
%Credit F. Fedderson
% returns the wavenumber of the gravity wave
% dispersion relation, by using newtons method
% the initial guess will be the shallow water wavenumber
ogf = f;
ogh = h;
omega = 2*pi*f;
g = 9.81;
k = omega./sqrt(g*h);
f = g*k.*tanh(k.*h) - omega.^2;
count = 0;
while max(abs(f))>1e-6
    dfdk = g*k.*h.*(sech(k.*h)).^2 + g*tanh(k.*h);
    k = k - f./dfdk;
    f = g*k.*tanh(k.*h) - omega.^2;
    if count > 1000
        break
    end
    count = count+1;
end
end
