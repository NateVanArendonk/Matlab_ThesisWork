function runups = editParkerRunupStack(runups)
%   runups = editParkerRunupStack(runups)
%
%  Allow user to edit a set of discrete runups to remove glitches.  When
%  editing, click using the left mouse button to add or improve bad points.
%   Any previous point that is within dtHalfWin will be replaced (this
%   value is displayed at the beginning of editing but is +/- 40% of the
%   average runups period.  

t = runups.tMax;                   % same here
r = runups.xMax;
ind = runups.maxInd;
T = runups.dn(:);
Nr = length(r);
dtHalfWin = runups.params.editHalfWinDT*(t(end)-t(1))/Nr;    % blank any within this of pick
foo = 'r';
if strcmp(foo,'r')
    disp('Left click new points, cr to finish')
    [x,y] = ginput;
    for i = 1: length(x)
        bad = find(abs(t-matlab2Epoch(x(i)))<dtHalfWin);
        if ~isempty(bad)    % eliminate replaced points
            t(bad) = nan;r(bad) = nan; ind(bad) = nan;
        end
        Nr = Nr+1;      % add new points
        dt = abs(T-x(i));
        t(Nr) = matlab2Epoch(x(i)); r(Nr)=y(i); ind(Nr) = find(dt == min(dt));
    end
    % now tidy up
    good = find(~isnan(t));
    t = t(good); r = r(good); ind = ind(good);
    % now sort
    [t,tInd] = sort(t);
    r = r(tInd); ind = ind(tInd);
    runups.tMax = t;       
    runups.xMax = r;
    runups.maxInd = ind;
    t = runups.tMax; r = runups.xMax;
    plot(epoch2Matlab(t),r,'go')

    runups = findRunupsExtremes(runups, runups.params.N);
    runups = findMaxForOverlapTime(runups, runups.RMaxFromBrightest.when);
    displayParkerResult(runups);
    
    dnRMax = T(runups.RMax.ind);
    dn2 = T(runups.R2.ind);
    dn5 = T(runups.R5.ind);
    dn10 = T(runups.R10.ind);
    dn50 = T(runups.R50.ind);
    plot(dnRMax, runups.RMax.x, 'co', 'markersize', 12, 'linewidth', 2)
    plot(dn2, runups.R2.x, 'm*', 'markersize', 12, 'linewidth', 2)
    plot(dn5, runups.R5.x, 'g*', 'markersize', 12, 'linewidth', 2)
    plot(dn10, runups.R10.x, 'go', 'markersize', 9, 'linewidth', 2)
    plot(dn50, runups.R50.x, 'go', 'markersize', 6, 'linewidth', 2)
    xb = runups.RMaxFromBrightest.x;
    t1 = runups.RMaxFromBrightest.when;
    t2 = t1+600;
    line([t1 t2], [xb xb], 'color', 'y')
    legend('runups', 'edited', 'RMax', 'R2', 'R5', 'brightest','R10','R50', 'location', 'southeast')
    runups.valid = 1;
else
    runups.valid = 0;
end

