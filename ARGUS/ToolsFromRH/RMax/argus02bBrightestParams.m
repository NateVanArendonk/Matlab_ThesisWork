% argus02b Brightest Params, 
% input for findBrightest analysis
% Specific support routines for a site
params.envInputProg = 'getDuckEnvirInfo';  % given epoch, return zt, Hs, fp
params.loadDEMProg = 'getLowTideDEM';      % given epoch, return DEM
params.findZEstProg = 'StockdonEquation';  % guess at RMax z-level

% Flag to use blue to red ratio (true) or grayshade (false)
params.useBRRatio = 1;       % 

% Search window params
params.halfWin = 10;         % half width of search window in meters
params.halfSubWin = 3;       % half width of final fit window in meters
params.halfSobelWidth = 0.5; % narrower width for initial sobel filtering
params.medFiltLength = 2;    % median filter length in meters
params.xRes = 0.05;          % interp into real space with this res.
params.guessLambda = 0.2;    % guess at e-folding length scale in meters

% sampling domain specs (extracting data from images)
params.ys = [0: 5: 1000]';    % sample transects
params.zs = 0;               % mapping level
params.minGood = 0.75;       % min fraction of pixels that must be on screen

% loessInterp specs for mapping RMaxDay output
params.loessInterp.xm = [60:1:140]; % map domain x (m)
params.loessInterp.ym = [0:5:1000]; % y domain (m)
params.loessInterp.lxy = [5 10];    % smoothing in x, y
params.loessInterp.sigGain = 2048; % sigma = sigGain/merit.  Huge = Kludge
params.loessInterp.cAxis = [-0.5 2.5];  % default colorscale
params.loessInterp.errcAxis = [-2 1];  % default log10(err) colorscale

% Miscellaneous
params.beta = 0.08;           % foreshore slope
params.betaSig = 0.005;       % initial error on slope guess
params.xsSig = 4;             % initial error on seed shoreline
params.maxInterpRMSE = 1;     % max err allowable in daily interp.
params.dyFitRange = 200;      % alongshore length of planar fit for new xs
params.QzParams = [0.05^2 0.05^2]; % process noise params for vertical
params.QxParams = [1^2 1^2]; % process noise params for xs
params.QBetaParams = [0.001^2 0.001^2]; % process noise params for beta

% Quality Control limits
params.QCLims = [50      200;   % x location
              -inf    inf;    % y location
              -inf    inf;    % z location
              0       100;    % camera number
              0.8     100;    % R2
              30      235;    % amplitude (actually range
              0.5      50;    % lambda, e-folding length scale
              0       225];   % I0, base brightness
    
% Debug options
params.debug.showShoreline = 1; % plot shoreline results on obliques
params.debug.showTransects = 0; % plot individual transects with pauses
params.debug.showMap = 1;       % show daily map of foreshore topo
params.debug.showAllKalman = 1; % show all Kalman components


