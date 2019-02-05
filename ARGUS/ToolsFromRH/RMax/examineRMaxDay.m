function examineRMaxDay(RMaxDay)
%   visualize output of RMax day analysis
%

% show individual collections
whichZStr = 'zEst';
showRMaxUpdateResults(RMaxDay, whichZStr)
showRMaxKalmanResults(RMaxDay, whichZStr)
showRMaxSeedAndBetaResults(RMaxDay)

