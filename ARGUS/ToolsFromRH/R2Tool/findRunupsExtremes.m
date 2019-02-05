function runups = findRunupsExtremes(runups, M)
%   runups = findRunupsExtremes(runups, M)
%
% take a set of runups (discrete maxima) and find the max as well as the R2
% and R5 values.  M is the number of columns in the stack image.

xMax = runups.xMax;
maxInd = runups.maxInd;
[~, ind] = sort(xMax);
runups.RMax.ind = maxInd(ind(1));
runups.RMax.x = xMax(ind(1));
if length(ind)>runups.params.minPeaksForStats      % don't find stats not enough peaks
    pick = ind(round(0.02*length(xMax)));
    runups.R2.ind = maxInd(pick);
    runups.R2.x = xMax(pick);
    pick = ind(round(0.05*length(xMax)));
    runups.R5.ind = maxInd(pick);
    runups.R5.x = xMax(pick);
    pick = ind(round(0.10*length(xMax)));
    runups.R10.ind = maxInd(pick);
    runups.R10.x = xMax(pick);
    pick = ind(round(0.5*length(xMax)));
    runups.R50.ind = maxInd(pick);
    runups.R50.x = xMax(pick);
else
    runups.R2.ind = M;
    runups.R2.x = nan;
    runups.R5.ind = M;
    runups.R5.x = nan;
    runups.R10.x = nan;
    runups.R10.x = nan;
    runups.R50.x = nan;
    runups.R50.x = nan;
end