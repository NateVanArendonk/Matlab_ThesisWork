function [wndspd,wnddir] = uv2compass(u,v)
%Converts from compass Wind Speed and Direction to U and V space
%   Note, this function works on meteorological data so it assumes that the
%   data is in meteorological convention.  This will output wind directions
%   as where the winds are coming from 

%   INPUTS: 
%       u: E/W component of winds as integers
%       v: N/S component of winds as integers 

% created 1/11/2019 - N.VanArendonk

%First calculate the windspeed 
wndspd = sqrt(u.^2 + v.^2);

%Now Calculate the wind direction 
wnddir = (180/pi).*atan2(v,u); % converts to a degree that is in polar coordinates 

% Convert to compass degrees 
wnddir = wrapTo360(90 - wnddir);

% Convert back to 'Coming From meteorological notation'
wnddir = wrapTo360(180 + wnddir);
end

