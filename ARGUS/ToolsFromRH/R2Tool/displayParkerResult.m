function     runups = displayParkerResult(runups)
%       runups = displayParkerResult(runupFullPathname);
%
% take a full pathname for a runup file from Parker's home, display the
% result including all of the extremes
if runups.valid
    load(FTPPath(runups.fn))
    dn = epoch2Matlab(T);
    [x, ind] = sort(XYZ(:,1));
    data = RAW(:,ind)';
    t = epoch2Matlab(runups.tMax);                   % same here
    r = runups.xMax;
    ind = runups.maxInd;
    x1 = mean(r)-5*std(r); x2 = mean(r)+3*std(r);
    good = find((x>=x1)&(x<=x2));
    figure(4); clf
    imagesc(dn,x(good),data(good,:))
    hold on; grid on; datetick('x')
    plot(t,r,'r*')
    % now plot the extremes
    plot(runups.dn(runups.RMax.ind), runups.RMax.x, 'co', 'markersize', 8, 'linewidth', 2)
    plot(runups.dn(runups.R2.ind), runups.R2.x, 'm*', 'markersize', 12, 'linewidth', 2)
    plot(runups.dn(runups.R5.ind), runups.R5.x, 'g*', 'markersize', 12, 'linewidth', 2)
    plot(epoch2Matlab(runups.RMax10Min.t), runups.RMax10Min.x, 'ko', 'markersize', 12, 'linewidth', 2)
    xb = runups.RMaxFromBrightest.x;
    t1 = epoch2Matlab(runups.RMaxFromBrightest.when);
    t2 = t1+600/(24*3600);
    legend('runups', 'RMax', 'R2', 'R5', 'brightest', 'location', 'southeast')
    line([t1 t2], [xb xb], 'color', 'b')
    title(datestr(epoch2Matlab(runups.when)))
end
