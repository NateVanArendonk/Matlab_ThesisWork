function [u,v] = compass2uv(wndspd,wnddir)
%Converts from compass Wind Speed and Direction to U and V space
%   Note, this function works on meteorological data so it assumes that the
%   data is in meteorological convention, i.e. the wind direction is in
%   reference to where the wind is coming from.  If it is not in this
%   convention, simply transform it with: NewDeg = wrapTo360(OldDeg + 180)

%   For this calculation angle will need to be converted from degrees to
%   radians 

%   INPUTS: 
%       wndspd: value or values corresponding to wind speeds (integers)
%       wnddir: value or values corresponding to wind direction (integers 
%       given as decimal degrees)

% created 1/11/2019 - N.VanArendonk

% Convert wind direction from radians to degrees 
wnddir = (pi/180)*wnddir;

% Convert from wind direction/speed to U and V space 
u = -wndspd.*sin(wnddir);
v = -wndspd.*cos(wnddir);
end

