function showRMaxUpdateResults(RMaxDay, whichZStr, figNum)
%   showRMaxUpdateResults(RMaxDay, whichZStr, {figNum})
%
% Plot update results from each collection (hourly or whatever) in an RMaxDay.
% Either 'zEst' or 'zDEM' can be chosen using whichZStr. The results will 
% be a 1 by 3 array of plots consisting of 1) the scatter plot, 2) the 
% loess-interped map, and 3) the rmse on the loess map.
% A non-default figure number can be entered.

if nargin < 3
    figNum = 10;
end

% extract all of the individual estimates into a large array
allxyzs = [];
for i = 1:length(RMaxDay.RMaxes)
    eval(['allxyzs = [allxyzs; RMaxDay.RMaxes(i).xy' whichZStr '];']);
end

% do a scatter3 plot of all the individual estimates
figure(figNum); clf; colormap(jet)
subplot(121)
scatter3(allxyzs(:,1), allxyzs(:,2), allxyzs(:,3), 4, allxyzs(:,3))
xm = RMaxDay.mapped.xm;
ym = RMaxDay.mapped.ym;
a = [min(xm) max(xm) min(ym) max(ym) RMaxDay.params.loessInterp.cAxis];
axis(a); view(2); caxis(a(5:6)); colorbar
hold on
plot(RMaxDay.seed.Kalman.xsKal, RMaxDay.params.ys, 'k')
dayStr = datestr(epoch2Matlab(RMaxDay.when), 2); 
xlabel('x (m)'); ylabel('y (m)'); title(['RMax scatter ' dayStr])

eval(['Z = RMaxDay.mapped.' whichZStr '.update.Z;']);
eval(['rmse = RMaxDay.mapped.' whichZStr '.update.rmse;']);
subplot(143)
imagesc(xm,ym,Z); caxis(RMaxDay.params.loessInterp.cAxis); colorbar
set(gca, 'ydir', 'nor')
%dn = epoch2Matlab(RMaxDay.when);
xlabel('x (m)'); ylabel('y (m)'); title([whichZStr ' loessUpdate (m)'])

subplot(144)
imagesc(xm,ym,log10(rmse));colorbar; caxis(RMaxDay.params.loessInterp.errcAxis)
set(gca, 'ydir', 'nor')
xlabel('x (m)'); ylabel('y (m)'); title('log10 rmse (m)')

