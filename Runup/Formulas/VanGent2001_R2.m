function R2 = VanGent2001_R2(Hs,To,beta)
% Calculates run up from VanGent, 2001: "Wave Run-Up on Dikes with Shallow Foreshores" - equation 4 for short waves: eq. 4
% Values for parameters c0 and c1 can be found in table 5 of paper. y is
% the reduction factor and is set to be 0.7.  

% Hs is Hsig in deep water
% To is peak period which will be converted to zeroth moment period 
% Beta is the slope 

% Constants
c0 = 1.55;
c1 = 5.4;
y = 0.7;


To = To/1.1; % Per guidance in VanderMeer.  Per VanGent "Using the peak wave period in predictions on wave run up 
% may lead to large inaccuracies for situations with shallow foreshores"
% (VanGent,2001)

% Calculate parameters needed for determining run up 
c2 = (0.25*(c1^2))/c0;
p = (0.5*c1)/c0;
Es = beta./(sqrt((2*pi/9.81)*(Hs./To.^2))); % Irribaren #

inds = Es <= p;
z2 = zeros(size(Hs));
z2(inds) = (y*Hs(inds)).*(c0*Es(inds));
z2(~inds)= (y*Hs(~inds)).*(c1-c2./Es(~inds));

R2 = z2;
end

