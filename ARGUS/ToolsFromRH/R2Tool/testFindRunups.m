% mess with runup extraction at Duck

clear
outPn = 'OUTPUT/';

pn = '/ftp/pub/argus02b/2015/cx/276_Oct.03/';
fns = dir([pn '*runup750*']);
runupParams

for i = 1: length(fns)
    fn = fns(i).name;
    load(FTPPath(fn))
    [x,ind] = sort(XYZ(:,1));
    data = double(RAW(:,ind));

    % start structure
    runups.when = T(1);
    runups.fn = fn;
    runups.params = params;

    % analyze and display
    [runups,I] = findRunups(x, T, data, runups);
    showRunupsOnStack(x, T, data, runups);
    if params.plot.showNoBackgroundVer
        showRunupsOnStackNoBackground(x,T,I,runups);
    end
    showRunupsHist(runups);
    pause

    % save
    p = parseFilename(fn);
    p.type = [p.type 'Proc'];
    fnOut = argusFilename(p);
    eval(['save ' outPn fnOut ' runups'])
end