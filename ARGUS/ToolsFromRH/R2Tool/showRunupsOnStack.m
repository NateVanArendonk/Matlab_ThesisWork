function showRunupsOnStack(x, t, data, runups)
%   showRunupsOnStack(x, t, data, runups)
%
%  Plot a runup stack with runup and picks

dn = runups.dn;
dnR = epoch2Matlab(runups.tMax);

NPanels = 3;
indPerPanel = ceil(length(dn)/NPanels);
figure(1); clf; colormap(gray)
for i = 1: NPanels
    subplot(1,NPanels,i); hold on
    inds = (i-1)*indPerPanel+1: min(i*indPerPanel, length(dn));
    imagesc(x,dn(inds),data(inds,:)); datetick('y')
    set(gca, 'ydir', 'rev'); xlabel('x (m)')
    xlim(runups.params.plot.xLims); 
    plot(runups.xr(inds),dn(inds))
    picks = find((dnR>dn(inds(1))) & (dnR<=dn(inds(end))));
    plot(runups.xMax(picks),dnR(picks),'r*')
    if ~isnan(runups.RMax.x)        % show extremes if they exist
        dnRMax = dn(runups.RMax.ind);
        dn2 = dn(runups.R2.ind);
        dn5 = dn(runups.R5.ind);
        plot(runups.RMax.x, dnRMax, 'co', 'markersize', 12, 'linewidth', 2)
        plot(runups.R2.x, dn2, 'm*', 'markersize', 12, 'linewidth', 2)
        plot(runups.R5.x, dn5, 'g*', 'markersize', 12, 'linewidth', 2)
    end
    if i==2
        if runups.valid
            title(safeString(runups.fn))
        else
            title([safeString(runups.fn) ' - NOT VALID'])
        end
    end
end
if ~isnan(runups.RMax.x)
    legend('runup', 'runups', 'RMax', 'R2', 'R5')
else
    legend('runup', 'runups')
end
