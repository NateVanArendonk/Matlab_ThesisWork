function showRunupsHist(runups)
%   showRunupsHist(runups)
%
%  plots a histogram and cumulative histogram of runups, but only if valid

figure(3); clf
if runups.valid
    xMax = runups.xMax;
    Nx = length(xMax);
    [N,xr] = hist(xMax);
    subplot(121)
    bar(xr, N/Nx);
    set(gca, 'xdir', 'rev')
    xlabel('xRunups (m)'); ylabel('Probability');
    p = parseFilename(runups.fn);
    title([p.station ', ' safeString(p.when) ', ' p.type])
    subplot(122)
    xMaxSort = sort(xMax, 1, 'descend');
    plot(xMaxSort, (1:Nx)/Nx)
    set(gca, 'xdir', 'rev'); grid on
    xlabel('xRunups (m)'); ylabel('cumulative prob.')
end
