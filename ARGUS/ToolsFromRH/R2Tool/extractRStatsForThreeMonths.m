% go through brightest output from four months in 2015, comparing against
% digitized runup

clear
DBConnect blossom.oce.orst.edu holman '' backup_argus

pn = '/home/ruby/users/holman/research/brightestWork/RMax/OUTPUT/';
outPn = 'OUTPUT/';
days = dir([pn '*RMax*']);

% for each day, find all the 850 values that are valid and do runup for
% them.
ymStr = '850';
iy = 171;       % you need to know this is the index for this ym
runupParams;

for day = 52: 52  %length(days)
    load([pn days(day).name]);
    for j = 1: length(RMaxDay.RMaxes)
        j
        if RMaxDay.RMaxes(j).results(iy,9) % if valid, get runup data
            runups = createEmptyRunupsStruct(params);
            runups.params = params;
            g = RMaxDay.RMaxes(j);
            runups.RMaxFromBrightest.x = g.results(iy,1);
            runups.RMaxFromBrightest.merit = g.results(iy,5)*g.results(iy,6);
            runups.RMaxFromBrightest.when = g.when;
            runups.ym = str2num(ymStr);
            runups.brightestParams = RMaxDay.params;
            r.when = g.when;
            r.ym = str2num(ymStr);
            fn = findArgusImages(r.when,'argus02b', ...
                'x',['runup' num2str(r.ym)], 'mat',100)
            if ~isempty(fn)
                runups = analyzeRunups(fn, runups);
                runups = findMaxForOverlapTime(runups, g.when);
                figure(4);
                plot(runups.RMax10Min.t, runups.RMax10Min.x,'ok', ...
                    'markersize', 12, 'linewidth', 2)
                pause(4)        % let user see result
                if runups.valid
                    p = parseFilename(fn);
                    p.type = [p.type 'Proc'];
                    fnOut = argusFilename(p);
                    eval(['save ' outPn fnOut ' runups'])
                end
            end
        end
    end
end

