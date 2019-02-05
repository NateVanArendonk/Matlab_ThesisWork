% create some figures for paper

clear
pn = '/home/ruby/users/holman/research/brightestWork/RMax/OUTPUT/';
fns = dir([pn '*Oct*RMaxDay*']);

i=25;       % a good choice
load([pn fns(i).name])
examineRMaxDay(RMaxDay)
subplot(143); set(gca, 'xtick', [60:40:140], 'xticklabel', [60:40:140])
subplot(144); set(gca, 'xtick', [60:40:140], 'xticklabel', [60:40:140])

allxyzs = [];
good = [];
whichZStr = 'zEst';
for i = 1:length(RMaxDay.RMaxes)
    eval(['allxyzs = [allxyzs; RMaxDay.RMaxes(i).xy' whichZStr '];']);
    good = [good; find(RMaxDay.RMaxes(i).results(:,9))];
end

% 3814 are good, max is 1.69, min is 0.13.

figure(13); clf
subplot(121)
imagesc(RMaxDay.mapped.xm, RMaxDay.mapped.ym, RMaxDay.mapped.zEst.Kalman.Z)
axis xy; xlabel('x (m)'); ylabel('y (m)'); colorbar
caxis([-1 3]); 
title(['Kalman topography, ' datestr(epoch2Matlab(RMaxDay.when),2)])
subplot(122)
imagesc(RMaxDay.mapped.xm, RMaxDay.mapped.ym, sqrt(RMaxDay.mapped.zEst.Kalman.P))
axis xy; xlabel('x (m)'); ylabel('y (m)'); colorbar
caxis([0 0.5]); 
title(['Kalman Error (m)'])