% go through Parkers runups one more time (my local copy), adding some
% things.
%   1.  Add reference to brightest image
%   2.  Collect environment stuff
%   3.  pick up shoreward max that were previously cut off by ylims
%   4.  Find other stats like R10, R20, R50

% modified to just fix the RMaxDay filenames
clear
pn = '/home/ruby/users/holman/research/RUNUP/discreteRunupTool/OUTPUTParker2/';
pnOut = '/home/ruby/users/holman/research/RUNUP/discreteRunupTool/OUTPUTParkerFinal2/';
fns = dir([pn '*mat']);
pnrm = '/home/ruby/users/holman/research/brightestWork/RMax/OUTPUT/';
fnsrm = dir([pnrm '*mat']);
foo = [fnsrm.name];foo = reshape(foo',[],length(fnsrm))';
trms = floor(epoch2Matlab(str2num(foo(:,1:10))));

for i = 1: length(fns)
    i
    load([pn fns(i).name]);
    ind1 = find(trms<=epoch2Matlab(runups.when),1,'last');
    load([pnrm fnsrm(ind1).name])
    tInDay = [RMaxDay.RMaxes.when];
    ind = find(tInDay >= runups.when, 1, 'first');
    indy = find(RMaxDay.RMaxes(ind).results(:,2) == runups.ym);
    runups.RMaxFromBrightest.RMaxDayfn = [pnrm fnsrm(ind1).name]; % assume cam 1
    runups.RMaxFromBrightest.env.zt = RMaxDay.RMaxes(ind).env.zt;
    runups.RMaxFromBrightest.env.Hs = RMaxDay.RMaxes(ind).env.Hs;
    runups.RMaxFromBrightest.env.fp = RMaxDay.RMaxes(ind).env.fp;
    runups.RMaxFromBrightest.env.beta = RMaxDay.RMaxes(ind).env.beta(indy);
    runups.RMaxFromBrightest.env.R2 = RMaxDay.RMaxes(ind).env.R2(indy);
    runups.RMaxFromBrightest.env.delz = RMaxDay.RMaxes(ind).env.delz(indy);
    runups.RMaxFromBrightest.env.delx = RMaxDay.RMaxes(ind).env.delx(indy);
%     runups = displayParkerResult(runups);
%     runups = editParkerRunupStack(runups);
    eval(['save ' pnOut fns(i).name ' runups'])
    %pause(5)
end

