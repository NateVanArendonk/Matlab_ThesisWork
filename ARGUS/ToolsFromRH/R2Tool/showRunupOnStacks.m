function showRunupOnStacks(x, t, data, runups)
%   showRunupOnStacks(x, t, data, runups)
%
%  Plot a runup stack with runup and picks

dn = epoch2Matlab(t);
dnR = epoch2Matlab(runups.tMax);
dnRMax = dnR(runups.RMax.ind);
dn2 = dnR(runups.R2.ind);
dn5 = dnR(runups.R5.ind);

NPanels = 3;
indPerPanel = ceil(length(dn)/NPanels);
figure(1); clf;
for i = 1: NPanels
    subplot(1,NPanels,i); hold on
    inds = (i-1)*indPerPanel+1: min(i*indPerPanel, length(dn));
    imagesc(x,dn(inds),data(inds,:)); datetick('y')
    set(gca, 'ydir', 'rev'); xlabel('x (m)')
    xlim(runups.params.plot.xLims); 
    plot(runups.xr(inds),dn(inds))
    picks = find((dnR>dn(inds(1))) & (dnR<=dn(inds(end))));
    plot(runups.xMax(picks),dnR(picks),'r*')
    plot(runups.RMax.x, dnRMax, 'wo')
    plot(runups.R2.x, dn2, 'y*')
    plot(runups.R5.x, dn5, 'g*')
    if i==2
        title(safeString(runups.fn))
    end
end
legend('runup', 'runups', 'RMax', 'R2', 'R5')

