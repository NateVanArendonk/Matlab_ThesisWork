% examineRunupsPicksOnStack
% load digitized runup and show digitized runups and stats

clear
pnIn = '/home/shusin-spare2/hamelp/runupWork/OUTPUT/';
pnSave = '/home/ruby/users/holman/research/RUNUP/discreteRunupTool/OUTPUTParker2/';
fns = dir([pnIn '*mat']);
addInds = [96:195];         % add the Sept days that Parker did on his last day

for i = addInds
    i
    figure(4); clf; colormap(gray)
    load([pnIn fns(i).name]);
    if runups.valid
        runups = displayParkerResult(runups);
        runups = editParkerRunupStack(runups);
        eval(['save ' pnSave fns(i).name ' runups'])
        pause(5)
    end
end

