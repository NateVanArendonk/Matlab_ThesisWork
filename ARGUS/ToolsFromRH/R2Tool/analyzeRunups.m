function runups = analyzeRunups(fn, runups)
%   runups = analyzeRunups(fn)
%
%  compute the individual runups and the swash for a runup, fn

outPn = 'OUTPUT/';
params = runups.params;

load(FTPPath(fn))
[x,ind] = sort(XYZ(:,1));
data = double(RAW(:,ind));

% start structure
runups.when = T(1);
runups.fn = fn;
runups.params = params;

% analyze and display
[runups,I] = findRunups(x, T, data, runups);
foo = 'r';
while(params.editAllowed & strcmp(foo, 'r') & (length(runups.tMax)>params.minNtMax))
    runups = editRunupStack(runups, x, T, I);
    if runups.valid
        foo = input('To re-do type r, otherwise hit cr to continue - ', 's');
    else
        return
    end
end
if params.plot.showRunupOnStacks
    showRunupsOnStack(x, T, data, runups);
end
if params.plot.showNoBackgroundVer
    showRunupsOnStackNoBackground(x,T,I,runups);
end
if params.plot.showHists
    showRunupsHist(runups);
end
    

% save
p = parseFilename(fn);
p.type = [p.type 'Proc'];
fnOut = argusFilename(p);
eval(['save ' outPn fnOut ' runups'])
