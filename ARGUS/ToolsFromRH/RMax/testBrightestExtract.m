% test brightest algorithms at Duck
% Development code.  
% Rob Holman, July 2018.

clear
addpath('/home/ruby/matlab/CILtest/CIL/argusDB2/');
DBConnect blossom.oce.orst.edu holman '' backup_argus

% pick a set of files (days) to digitize, create a seed structure
mon = 'Nov';
pn = '/ftp/pub/argus02b/2015/c1/';
days = dir([pn '*' mon '*']);
fns = dir([pn days(1).name '/*bright*']);
fnPick = fns(round(length(fns)/2)).name;        % pick a mid-day example seed
n = parseFilename(fnPick);                      % for later use
camInfo = findCamInfo(n.station, str2num(n.time));
camNums = [camInfo.num];
RMaxDay = createFirstRMaxDayRecord(fnPick, camInfo);        % create a seed day

% cycle through the days of the month
tic
for dayNum = 1: length(days)
    fns = dir([pn days(dayNum).name '/*bright*']);  % get a day list of C1 fns
    if length(fns)>=1                           % analyze if any images
        for i = 1: length(fns)                      % cycle through the images
            fnList(1,:) = fns(i).name;
            n = parseFilename(fns(i).name);
            for c = 2: length(camNums)
                n.camera = c;
                foo = findArgusImages(n, 60);  % restrict to near synchronous
                if size(foo,1) == 1
                    fnList(c,:) = foo;
                else
                    fnList(c,:) = repmat(' ', size(fnList(1,:)));
                end
            end
            % compute the results for that image. Use B/R ratio (last arg)
            RMaxes = findRMax(fnList, RMaxDay, camInfo, RMaxDay.params.useBRRatio);    
            RMaxes.xyzEst = findxyzEst(RMaxes, camInfo);
            % now load a DEM and find xyzDEM
            eval(['d = ' RMaxDay.params.loadDEMProg '(RMaxes.when);']);
            if ~isempty(d.fullPn)
                RMaxes.xyzDEM = findxyzDEM(RMaxes, camInfo, d);
                RMaxes.DEMpn = d.fullPn;
            else
                RMaxes.xyzDEM = nan(size(RMaxes.xyzEst));
                RMaxes.DEMpn = [];
            end
            RMaxDay.RMaxes(i) = RMaxes;
        end
        RMaxDay.RMaxes = RMaxDay.RMaxes(1:length(fns));  % remove left overs
    end
    RMaxDay = makeMapsRMaxDay(RMaxDay);      % assemble the full day and save
    p.time = RMaxDay.when;
    p.station = n.station;
    p.camera = 'x';
    p.type = 'RMaxDay';
    p.format = 'mat';
    fnOut = argusFilename(p);
    eval(['save OUTPUT/' fnOut ' RMaxDay'])
    examineRMaxDay(RMaxDay)                 % examine if useful
end
Tau = toc
computeTimePerDay = Tau/length(days)
