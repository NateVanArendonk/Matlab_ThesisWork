function RMaxDay = makeMapsRMaxDay(RMaxDay)
%   RMaxDay = makeMapsRMaxDay(RMaxDay)
%
% Given a full day of RMax estimates, 
% compute 
%   map     - a 2D loess-interped topography
%   seed    - a shoreline vector at z=0 for seeding the search
%   beta    - a beach slope vector for estimated R2.
% Each quantity is returned as an update (particular to that collection)
% and a Kalman (smoothed with a Kalman filter) version.  Support data are
% contained with each structure.  

params = RMaxDay.params;
RMax = RMaxDay.RMaxes;
allxyzs = [];
merit = [];
etList = [RMaxDay.RMaxes.when];
RMaxDay.when = mean(etList);
RMaxDay.whenPrevious = RMaxDay.mapped.zEst.Kalman.when;
dt = (RMaxDay.when - RMaxDay.whenPrevious) / (24*3600); 

%% 1. Collect data from the full day.  Create an ad hoc merit as R-sq * ht.
for i = 1: length(RMax)
    r = RMax(i).results;
    good = find(r(:,9));
    allxyzs = [allxyzs; RMax(i).xyzEst(good,:)];
    merit = [merit; [r(good,5).*r(good,6)]];
    Hs(i) = RMax(i).env.Hs;
end
% create a bogus CI from merit (product of R2 and amp)
sigZ = params.loessInterp.sigGain./merit;

%% 2.  Do a loessInterped daily map in update and Kalman format
xm = params.loessInterp.xm;
ym = params.loessInterp.ym;
[X,Y] = meshgrid(xm, ym);       % design grid
warning off
[Z, rmse, ~, ~, sk,dof] = loessInterp([allxyzs(:,1) allxyzs(:,2)], ...
    allxyzs(:,3), sigZ, [X(:) Y(:)], params.loessInterp.lxy);
warning on
Z = reshape(Z, size(X));
rmse = reshape(rmse, size(X));
bad = find(rmse>params.maxInterpRMSE);         % remove bad data
Z(bad) = nan;
RMaxDay.mapped.xm = xm;            % save interp-ed map
RMaxDay.mapped.ym = ym;
RMaxDay.mapped.zEst.update.Z = Z;
RMaxDay.mapped.zEst.update.rmse = rmse; % square this for Kalman R

% now Kalman smooth the interp-ed estimates
HBar = mean(Hs);            % find day average Hs
POld = RMaxDay.mapped.zEst.Kalman.P;
ZOld = RMaxDay.mapped.zEst.Kalman.Z;    % previous day's map
Q = findQ(RMaxDay.params.QzParams, HBar, dt);  % find process error (depends on Hs)
[ZNew, PNew, K] = KalmanFilter(Z, ZOld, POld, rmse.^2, Q);
RMaxDay.mapped.zEst.Kalman.Z = ZNew;    % update to new map
RMaxDay.mapped.zEst.Kalman.ZOld = ZOld;
RMaxDay.mapped.zEst.Kalman.P = PNew;
RMaxDay.mapped.zEst.Kalman.POld = POld;
RMaxDay.mapped.zEst.Kalman.Q = Q;
RMaxDay.mapped.zEst.Kalman.K = K;
RMaxDay.mapped.zEst.Kalman.when = RMaxDay.when;

%% 3. Now do the same entire process using the DEM data.
for i = 1: length(RMax)
    r = RMax(i).results;
    good = find(r(:,9));
    allxyzs = [allxyzs; RMax(i).xyzDEM(good,:)];
    merit = [merit; [r(good,5).*r(good,6)]];
    Hs(i) = RMax(i).env.Hs;
end
% create a bogus CI from merit (product of R2 and amp)
sigZ = params.loessInterp.sigGain./merit;

%  Do a loessInterped daily map in update and Kalman format
warning off
[Z, rmse, ~, ~, sk,dof] = loessInterp([allxyzs(:,1) allxyzs(:,2)], ...
    allxyzs(:,3), sigZ, [X(:) Y(:)], params.loessInterp.lxy);
warning on
Z = reshape(Z, size(X));
rmse = reshape(rmse, size(X));
bad = find(rmse>params.maxInterpRMSE);         % remove bad data
Z(bad) = nan;
RMaxDay.mapped.zDEM.update.Z = Z;
RMaxDay.mapped.zDEM.update.rmse = rmse; % square this for Kalman R

% now Kalman smooth the interp-ed estimates
POld = RMaxDay.mapped.zDEM.Kalman.P;
ZOld = RMaxDay.mapped.zDEM.Kalman.Z;    % previous day's map
[ZNew, PNew, K] = KalmanFilter(Z, ZOld, POld, rmse.^2, Q);
RMaxDay.mapped.zDEM.Kalman.Z = ZNew;    % update to new map
RMaxDay.mapped.zDEM.Kalman.ZOld = ZOld;
RMaxDay.mapped.zDEM.Kalman.P = PNew;
RMaxDay.mapped.zDEM.Kalman.POld = POld;
RMaxDay.mapped.zDEM.Kalman.Q = Q;
RMaxDay.mapped.zDEM.Kalman.K = K;
RMaxDay.mapped.zDEM.Kalman.when = RMaxDay.when;

%% 4.  Now find a new shoreline seed and beach slope for the next day.  
% Do planar fits over a longshore range of params.dyFitRange for each ym. 
% Do this on the zEst Kalman data.
dy = median(diff(params.ys));
dyFitInd = round(params.dyFitRange/2/dy);
R = nan(size(params.ys));
xOld = RMaxDay.seed.Kalman.xsKal;  % grab previous shoreline
for i = 1:length(ym)
    inds = [max(1,i-dyFitInd): min(length(ym),i+dyFitInd)];
    x = X(inds,:); y = Y(inds,:); z = ZNew(inds,:);
    good = find(~isnan(z));
    beta0 = [params.beta 0 min(x(good))];
    if length(good)>20
        xy = [x(good) y(good)-ym(i)];   % center on transect of interest
        [betaAll, Resid, J, covB] = nlinfit(xy, z(good), @TwoDPlaneFit, beta0);
        zs = TwoDPlaneFit(betaAll,xy);
        xs(i) = -betaAll(3)/betaAll(1);
        beta(i) = betaAll(1);
        dz = zs - z(good);
        rmse(i) = sqrt(mean(dz.^2));
        CI = nlparci(betaAll, Resid, 'Jacobian', J);
        betaSig(i) = (CI(1,2)-CI(1,1))/2;
        R(i) = ((CI(3,2)-CI(3,1))/2/betaAll(1)).^2;    % approx error on xs
    else
        xs(i) = xOld(i);    % fake bad data by maintaining old value
        rmse(i) = nan;
        R(i) = 10000;        % use large R
        beta(i) = params.beta;  % default
        betaSig(i) = nan;
    end
end
RMaxDay.seed.update.xsAtz0 = xs(:);
RMaxDay.seed.update.rmse = sqrt(R(:));

% Kalman filter shoreline, xs, saving in Kalman.xsKal
dt = RMaxDay.when - RMaxDay.seed.Kalman.when;     % how long since previous estimate
dt = abs(dt)/(3600*24);     % convert to days
Q = findQ(params.QxParams, HBar, dt);
POld = RMaxDay.seed.Kalman.PNew;        % previous new result
R = R(:);
xs = xs(:);
[xsAtz0, P, K] = KalmanFilter(xs, xOld, POld, R, Q);
RMaxDay.seed.Kalman.xsKal = xsAtz0;       % update xs
RMaxDay.seed.Kalman.xsOld = xOld;
RMaxDay.seed.Kalman.PNew = P;     % update P
RMaxDay.seed.Kalman.POld = POld;     
RMaxDay.seed.Kalman.Q = Q;
RMaxDay.seed.Kalman.K = K;
RMaxDay.seed.Kalman.when = RMaxDay.when;

% beach slope results including Kalman.
RMaxDay.beta.update.beta = beta(:);
RMaxDay.beta.update.betaSig = betaSig(:);
betaOld = RMaxDay.beta.Kalman.betaKal;
POld = RMaxDay.beta.Kalman.PNew;        % previous
Q = repmat(findQ(params.QBetaParams,HBar,dt), length(ym),1);
R = betaSig.^2; R = R(:);
beta = beta(:);
[betaNew, P, K] = KalmanFilter(beta, betaOld, POld, R, Q);
RMaxDay.beta.Kalman.betaKal = betaNew;
RMaxDay.beta.Kalman.betaOld = betaOld;
RMaxDay.beta.Kalman.PNew = P;
RMaxDay.beta.Kalman.POld = POld;
RMaxDay.beta.Kalman.Q = Q;
RMaxDay.beta.Kalman.K = K;
RMaxDay.beta.Kalman.when = RMaxDay.when;


function z = TwoDPlaneFit(beta, xy)
%
%   2D planar fit z = -(beta(1)*x + beta(2)*y + beta(3)).  Note z = -h.
%
    z = -(beta(1)*xy(:,1) + beta(2)*xy(:,2) + beta(3));
    
