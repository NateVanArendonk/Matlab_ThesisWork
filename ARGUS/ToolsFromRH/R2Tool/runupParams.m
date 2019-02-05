% peak runup params

params.NFiltDark = 3;           % fast approach for darkening
params.NFiltLight = 30;         % slow approach to brightening
params.minIMax = 40;            % min mean stack brightness above background
params.upRushFraction = 0.4;    % min edge intensity for uprush
params.downRushFraction = 0.2;  % min edge intensity for downrush
params.minDt = 1.4;             % min dt in seconds for spike definition
params.maxDropThruFract = 0.01; % allow small number of drop throughs
params.minPeaksForStats = 50;   % needed runups for stats to be useful
params.editHalfWinDT = 0.4;     % fraction of mean runups period to edit
params.N = 2048;                % number of points in time series
params.minNtMax = 30;           % min number of maxima to continue

% plotting params
params.plot.xLims = [70 130];
params.plot.showRunupOnStacks = 0;
params.plot.showNoBackgroundVer = 0;    % show stack with IBack removed
params.plot.showHists = 0;
params.editAllowed = 1;

