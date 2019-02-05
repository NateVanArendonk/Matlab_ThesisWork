function r = createEmptyRunupsStruct(params)
%   runups = createEmptyRunupsStruct(params);
%
%  Self evident

% common material
N = params.N;        
r.when = [];
r.fn = [];
r.ym = nan;
r.params = [];
r.brightestParams = nan;
r.valid = 0;

% runup time series - length of stack
r.dn = nan(N,1);
r.xr = nan(N,1);
r.rInd = nan(N,1);

% runups - discrete maxima from runup time series
r.tMax = nan;
r.xMax = nan;
r.maxInd = nan;
r.IMax = nan;
% runups - statistical extremes of runups
r.RMax.ind = nan;
r.RMax.x = nan;
r.R2.ind = nan;
r.R2.x = nan;
r.R5.ind = nan;
r.R5.x = nan;

% results from or corresponding to brightest images
r.RMaxFromBrightest.when = nan;
r.RMaxFromBrightest.x = nan; % extracted from brightest image
r.RMaxFromBrightest.merit = nan;
r.RMax10Min.ind = nan;      % results from stack for brightest sample time
r.RMax10Min.x = nan;
r.RMax10Min.t = nan;

