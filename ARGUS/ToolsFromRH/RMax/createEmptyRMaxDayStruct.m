function RMaxDay = createEmptyRMaxDayStruct(NCams, params)
%   RMaxDay = createEmptyRMaxDayStruct(NCams, params)
%
%  Create the basic RMaxDay structure.  You need to specify the number of
%  cameras, NCams, and the params structure (which contains info about the
%  size of various fields).

% basic elements
RMaxDay.when = nan;
RMaxDay.whenPrevious = nan;
RMaxDay.version = getRMaxVersion;
RMaxDay.params = params;
RMaxDay.camNums = nan(NCams,1);
RMaxDay.xyzCams = nan(NCams,3);
Ny = length(params.ys);
nan1D = nan(Ny,1);

% now start creating substructures.
env.et = nan;
env.zt = nan;
env.Hs = nan;
env.fp = nan;
env.REst = nan(Ny,1);
env.delz = nan(Ny,1);
env.delx = nan(Ny,1);

% RMaxes, the results for each image set
RMaxes.when = nan;
RMaxes.whenDone = nan;
RMaxes.env = env;
RMaxes.fnList = [];
RMaxes.results = nan(Ny, 9);
RMaxes.xyzEst = nan(Ny,3);
RMaxes.xyzDEM = nan(Ny,3);
RMaxes.DEMpn = [];

% mapped, the interpolated map results
mapped.xm = params.loessInterp.xm;
mapped.ym = params.loessInterp.ym;
Nxm = length(mapped.xm); Nym = length(mapped.ym);
nan2DMap = nan(Nym, Nxm);
mapped.zEst.update.Z = nan2DMap;        % map results from xyzEst
mapped.zEst.update.rmse = nan2DMap;
mapped.zEst.Kalman.Z = nan2DMap;
mapped.zEst.Kalman.P = nan2DMap;
mapped.zEst.Kalman.ZOld = nan2DMap;
mapped.zEst.Kalman.POld = nan2DMap;
mapped.zEst.Kalman.Q = nan2DMap;
mapped.zEst.Kalman.K = nan2DMap;
mapped.zEst.Kalman.when = nan;
mapped.zDEM.update.Z = nan2DMap;         % map results from DEM
mapped.zDEM.update.rmse = nan2DMap;
mapped.zDEM.Kalman.Z = nan2DMap;
mapped.zDEM.Kalman.P = nan2DMap;
mapped.zDEM.Kalman.ZOld = nan2DMap;
mapped.zDEM.Kalman.POld = nan2DMap;
mapped.zDEM.Kalman.Q = nan2DMap;
mapped.zDEM.Kalman.K = nan2DMap;
mapped.zDEM.Kalman.when = nan;

% shoreline seed. Map using zEst Kalman
seed.update.xsAtz0 = nan1D;
seed.update.rmse = nan1D;
seed.Kalman.xsOld = nan1D;
seed.Kalman.xsKal = nan1D;
seed.Kalman.POld = nan1D;
seed.Kalman.PNew = nan1D;
seed.Kalman.Q = nan1D;
seed.Kalman.K = nan1D;
seed.Kalman.when = nan;

% beach slope results.  Base on zEst Kalman
beta.update.beta = nan1D;
beta.update.betaSig = nan1D;
beta.Kalman.betaOld = nan1D;
beta.Kalman.betaKal = nan1D;
beta.Kalman.POld = nan1D;
beta.Kalman.PNew = nan1D;
beta.Kalman.Q = nan1D;
beta.Kalman.K = nan1D;
beta.Kalman.when = nan;

% now assemble all into full RMaxDay structure
RMaxDay.RMaxes = RMaxes;        % note that this is just one element
RMaxDay.mapped = mapped;
RMaxDay.seed = seed;
RMaxDay.beta = beta;



