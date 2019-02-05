function showRMaxSeedAndBetaResults(RMaxDay)
%   showRMaxSeedAndBetaResults(RMaxDay)
%
% plots the seed, xsAtz0, and the beach slope, beta, for an RMaxDay

ys = RMaxDay.params.ys;

% plot the shoreline seed update, Kalman and Kalman gain
figure(12); clf
subplot(131); hold on
xs = RMaxDay.seed.update.xsAtz0;
ci = RMaxDay.seed.update.rmse;
plot(xs, ys, 'k')
hold on
plot(xs+ci, ys, 'g')
plot(xs-ci, ys, 'g')
xlabel('x (m)'); ylabel('y (m)'); 
dn = epoch2Matlab(RMaxDay.when);
title(['xs update ' datestr(dn,2)]); grid on

subplot(132); hold on
xsKal = RMaxDay.seed.Kalman.xsKal;
ciKal = sqrt(RMaxDay.seed.Kalman.PNew);
plot(xsKal, ys, 'k')
hold on
plot(xsKal+ciKal, ys, 'g')
plot(xsKal-ciKal, ys, 'g')
xlabel('x (m)'); ylabel('y (m)'); 
dn = epoch2Matlab(RMaxDay.when);
title(['xs Kalman ' datestr(dn,2)]); grid on

subplot(133)
K = RMaxDay.seed.Kalman.K;
plot(K, ys, 'k')
title('Kalman Gain, K'); grid on

figure(13); clf
subplot(131)
beta = RMaxDay.beta.update.beta;
ci = RMaxDay.beta.update.betaSig;
plot(beta, ys, 'k')
hold on
plot(beta+ci, ys, 'g')
plot(beta-ci, ys, 'g')
xlabel('x (m)'); ylabel('y (m)'); 
title('beach slope update'); grid on

subplot(132)
betaKal = RMaxDay.beta.Kalman.betaKal;
betaci = sqrt(RMaxDay.beta.Kalman.PNew);
plot(betaKal, ys, 'k')
hold on
plot(betaKal+betaci, ys, 'g')
plot(betaKal-betaci, ys, 'g')
xlabel('x (m)'); ylabel('y (m)'); 
title('beach slope Kalman'); grid on

subplot(133)
K = RMaxDay.beta.Kalman.K;
plot(K, ys, 'k')
title('Kalman Gain, K'); grid on
