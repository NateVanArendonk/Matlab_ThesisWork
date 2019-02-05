function examineDetailsOfzEstKalman(RMaxDay)
%   examineDetailsOfzEstKalman(RMaxDay)
%
% plot all the components of a Kalman filter process for zEst for any
% RMaxDay

xm = RMaxDay.mapped.xm; ym = RMaxDay.mapped.ym;
figure(20); clf
subplot(151)
imagesc(xm, ym, RMaxDay.mapped.zEst.update.Z)
title('new Z'); colorbar; caxis([-0.5 2.5]); axis xy
subplot(152)
imagesc(xm, ym, RMaxDay.mapped.zEst.Kalman.ZOld)
title('old Kalman Z'); colorbar; caxis([-0.5 2.5]); axis xy
subplot(153)
imagesc(xm, ym, RMaxDay.mapped.zEst.Kalman.Z)
title('new Kalman Z'); colorbar; caxis([-0.5 2.5]); axis xy
subplot(154)
imagesc(xm, ym, log10(RMaxDay.mapped.zEst.update.rmse))
title('log10 new rmse'); colorbar; caxis([-2 1]); axis xy
subplot(155)
imagesc(xm, ym, log10(sqrt(RMaxDay.mapped.zEst.Kalman.POld)))
Q = RMaxDay.mapped.zEst.Kalman.Q;
QStr = num2str(sqrt(Q), '%.3f');
title(['log10 old rmse, Q = ' QStr]); colorbar; caxis([-2 1]); axis xy

figure(21); clf
subplot(141)
imagesc(xm, ym, log10(RMaxDay.mapped.zEst.update.rmse))
title('log10 new rmse'); colorbar; caxis([-2 1]); axis xy
subplot(142)
imagesc(xm, ym, log10(sqrt(RMaxDay.mapped.zEst.Kalman.POld+Q)))
title('log10 aPriori rmse'); colorbar; caxis([-2 1]); axis xy
subplot(143)
imagesc(xm, ym, RMaxDay.mapped.zEst.Kalman.K)
title('K'); colorbar; caxis([0 1]); axis xy
subplot(144)
imagesc(xm, ym, log10(sqrt(RMaxDay.mapped.zEst.Kalman.P)))
title('log10 aPosteriori rmse'); colorbar; caxis([-2 1]); axis xy
