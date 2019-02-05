function [zOut, PNew, K] = KalmanFilter(zNew, zOld, P, R, Q)
%   [zOut, P, K] = KalmanFilter(zNew, zOld, P, R, Q)
%
% Kalman update an old estimate zOld with variance P using a new
% measurement zNew with variance R, allowing for process error Q.  z can be
% any number of dimensions since the calculations are done point by point.
% The Kalman variance, P, is updated and returned as well as the Kalman
% gain.

P = P + Q;              % add process error update
K = P ./ (P+R);         % find Kalman gains
zOut  = zOld + K .* (zNew-zOld);
PNew = (1 - K) .* P;       % adjust variance.

% take care of nan's
badNew = isnan(zNew);
badOld = isnan(zOld);

badNewInd = find(badNew);
zOut(badNewInd) = zOld(badNewInd);
PNew(badNewInd) = P(badNewInd);
K(badNewInd) = 0;          % forced K to be consistent with nans

badOldInd = find(badOld);
zOut(badOldInd) = zNew(badOldInd);
PNew(badOldInd) = R(badOldInd);
goodNew = find(~badNew & badOld);
K(goodNew) = 1;          % same