function runups = editRunupStack(runups,x,T,I)
%   runups = editRunupStack(runups,x,T,I)
%
%  Allow user to edit a set of discrete runups to remove glitches.  When
%  editing, click using the left mouse button to add or improve bad points.
%   Any previous point that is within dtHalfWin will be replaced (this
%   value is displayed at the beginning of editing but is +/- 40% of the
%   average runups period.  

figure(4); clf
imagesc(T,x, I'); colormap(gray)   % make time local
hold on; grid on
%plot(T, runups.xr)
t = runups.tMax;                   % same here
r = runups.xMax;
ind = runups.maxInd;
plot(t,r,'r*')
xBar = mean(r); xSig = std(r); nSig = 3;
xMin = xBar-nSig*xSig; xMax = xBar+nSig*xSig;
ylim([xMin xMax])

% toss anything outside of nSig standard deviations
good = find((r>xMin) & (r<xMax));
r = r(good); ind = ind(good); t = t(good);
Nr = length(r);
dtHalfWin = runups.params.editHalfWinDT*(t(end)-t(1))/Nr    % blank any within this of pick
a = axis;           % add a scale for +/- dtHalfWin
x0 = a(1)+0.05*(a(2)-a(1)); y0 = a(3)+0.1*(a(4)-a(3));
line([x0-dtHalfWin x0+dtHalfWin],[y0 y0], 'color', 'w');
line([x0 x0], [y0-0.02*(a(4)-a(3)), y0+0.02*(a(4)-a(3))], 'color', 'w')
foo = input('Type <CR> to edit, n to dump this file - ','s');
if ~strcmp(foo,'n')
    disp('Left click new points, cr to finish')
    [x,y] = ginput;
    for i = 1: length(x)
        bad = find(abs(t-x(i))<dtHalfWin);
        if ~isempty(bad)    % eliminate replaced points
            NBad = length(bad);
            t(bad) = nan;r(bad) = nan; ind(bad) = nan;
        end
        Nr = Nr+1;      % add new points
        dt = abs(T-x(i));
        t(Nr) = x(i); r(Nr)=y(i); ind(Nr) = find(dt == min(dt));
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
    plot(t,r,'go')

    runups = findRunupsExtremes(runups, size(I,2));
    edn = matlab2Epoch(runups.dn(:));
    dnRMax = edn(runups.RMax.ind);
    dn2 = edn(runups.R2.ind);
    dn5 = edn(runups.R5.ind);
    plot(dnRMax, runups.RMax.x, 'co', 'markersize', 12, 'linewidth', 2)
    plot(dn2, runups.R2.x, 'm*', 'markersize', 12, 'linewidth', 2)
    plot(dn5, runups.R5.x, 'g*', 'markersize', 12, 'linewidth', 2)
    xb = runups.RMaxFromBrightest.x;
    t1 = runups.RMaxFromBrightest.when;
    t2 = t1+600;
    line([t1 t2], [xb xb], 'color', 'y')
    legend('runups', 'junk', 'junk', 'edited', 'RMax', 'R2', 'R5', 'brightest')
    runups.valid = 1;
else
    runups.valid = 0;
end

