function [runups, I] = findRunups(x, t, data, runups)
%   [runups, I] = findRunups(x, epochTime, data, runups)
%
%  find the runup time series from the stack, data(t,x), using params
%  specified in runups.params.  x in meters, t in epoch seconds, data is
%  stack as matlab doubles, gray shade.  The stack with the back removed is
%  also returned as an optional output, I.

params = runups.params;
I = removeBackground(data, params);
M = size(I,2);
[xr, rInd, IMax] = findRunupEngine(x, I, params);
% now find the maxima
[tMax, xMax, maxInd] = findRunupMaxima(t, xr, params);

% check for drop throughs
bad = find(rInd == size(I,2));
if (length(bad) > params.maxDropThruFract*size(I,1))   
    runups.valid = 0;
else
    runups.valid = 1;
end

runups.dn = epoch2Matlab(t);
runups.xr = xr(:);
runups.rInd = rInd(:);
runups.tMax = tMax(:);
runups.xMax = xMax(:);
runups.maxInd = maxInd(:);
runups.IMax = IMax;

runups = findRunupsExtremes(runups,M);

