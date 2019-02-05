% compare Jenna's runup data with mine.

clear
Jpn = 'JennaData/077_Mar.18_timeseries/';
fns = dir([Jpn '*r750*']);
for i = 1: length(fns)
    fn = fns(i).name;
    load([Jpn fn])

    dn = epoch2Matlab(runup.epoch);
    rJx = runup.xyz(:,1);

    runupParams
    p = parseFilename(fn);
    p.format = 'mat';
    p.type = ['runup' p.type(2:end)];
    sn = findArgusImages(p);

    runups = analyzeRunups(sn);

    figure(4);clf; colormap(gray)
    plot(dn, rJx); hold on
    plot(runups.dn, runups.xr, 'r--')
    ylim(params.plot.xLims)
    datetick('x')

    pause
end