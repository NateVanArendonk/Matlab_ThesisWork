function showRMaxKalmanResults(RMaxDay, whichZStr, figNum)
%   showRMaxKalmanResults(RMaxDay, whichZStr, {figNum})
%
%  Plot Kalman foreshore topography map along with error bars (log10
%  scaled) for any RMaxDay structure.  You can specify a figure number if
%  you don't like the default.  WhichZStr allows you to chose the vertical
%  reference, either 'zEst', or 'zDEM'.

if nargin<3
    figNum = 11;
end

xm = RMaxDay.mapped.xm;
ym = RMaxDay.mapped.ym;
eval(['ZUp = RMaxDay.mapped.' whichZStr '.update.Z;']);
eval(['ZKal = RMaxDay.mapped.' whichZStr '.Kalman.Z;']);
eval(['rmse = sqrt(RMaxDay.mapped.' whichZStr '.Kalman.P);']);
eval(['K = RMaxDay.mapped.' whichZStr '.Kalman.K;']);

figure(figNum); clf; colormap(jet)
subplot(141)
imagesc(xm,ym,ZUp); caxis(RMaxDay.params.loessInterp.cAxis); colorbar
set(gca, 'ydir', 'nor')
dn = epoch2Matlab(RMaxDay.when);
xlabel('x (m)'); ylabel('y (m)'); title([whichZStr ' (m),  ' datestr(dn,2), ', update'])

subplot(142)
imagesc(xm,ym,ZKal); caxis(RMaxDay.params.loessInterp.cAxis); colorbar
set(gca, 'ydir', 'nor')
dn = epoch2Matlab(RMaxDay.when);
xlabel('x (m)'); ylabel('y (m)'); title([whichZStr ' (m), ' datestr(dn,2), ', Kalman'])

subplot(143)
imagesc(xm,ym,log10(rmse));colorbar; caxis(RMaxDay.params.loessInterp.errcAxis)
set(gca, 'ydir', 'nor')
xlabel('x (m)'); ylabel('y (m)'); title('log10 rmse (m)')

subplot(144)
imagesc(xm,ym,K); colorbar; caxis([0 1])
set(gca, 'ydir', 'nor')
xlabel('x (m)'); ylabel('y (m)'); title('K')

