
function lims = argus02bQCLimits
% QCLimits for brightest results from argus02b
%  criteria for good brightest results for each variable

lims = [50      200;    % x location
        -inf    inf;    % y location
        -inf    inf;    % z location
        0       100;    % camera number
        0.8     100;    % R2
        25      255;    % amplitude
        0.5     150;     % lambda, e-folding length scale
        0       200];   % I0, base brightness
