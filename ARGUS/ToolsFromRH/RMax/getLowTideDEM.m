function d = getLowTideDEM(e)
%   d = getLowTideDEM(epochTime)
%
%   get the Lidar survey from Duck that corresponds to the lowest tide in a
%   +/- 12 hour period centered at epochTime.

twelveHour = 12*3600;
[et, zt] = DBGetPredictedTide('argus02b', [e-twelveHour e+twelveHour]);
[~, ind] = min(zt);
eLow = et(ind);
d = getClosestDEM(eLow);
