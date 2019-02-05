function angleOut = smallestSignedAngleBetween(angle1,angle2)
%This function returns the smallest angle between to angles, keeping the
%sign of the angle, necessary for computing alongshore radiative stress of waves   

% For Example. angle1 = 148
%              angle2 = 300 
%              angleOut is -152

% Inputs:
%   angle1: First Angle in Degrees
%   angle2: Second Angle in Degrees 

% Outputs:
%   angleOut: Smallest angle between the two, with the correct sign 


% Notes:  If two angles are X and Y then one of the angles between them is
% equal to the abs(x -y).  The Other angle is (2*pi)-abs(x-y)
% Source: https://stackoverflow.com/questions/1878907/the-smallest-difference-between-2-angles
% God Bless you Stack overflow 

% Function modified from above link by N. VanArendonk - 1/18/2019

Tau = 2*pi; 

d1 = deg2rad(angle1); % Convert first angle to radians
d2 = deg2rad(angle2); % Convert second angle to radians 

a = rad2deg(mod((d1 - d2),Tau)); % Find first angle difference within 0 and 2pi
b = rad2deg(mod((d2 - d1),Tau)); % Find second angle difference within 0 and 2pi

if a < b % If the second angle is smaller 
    angleOut = -a; % Return that and flip the sign  
else
    angleOut = b; % Otherwise return the other angle 
end

end

