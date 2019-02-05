% examine results for which there are large brightest error

clear
pnParkerOutFinal2 = '/home/ruby/users/holman/research/RUNUP/discreteRunupTool/OUTPUTParkerFinal2/';


fns = dir([pnParkerOutFinal2 '*runup*'])
ShusinWC = length(fns)
foo = [fns.name];
foo = reshape(foo',[],ShusinWC)';
shUnique = unique(floor(epoch2Matlab(str2num(foo(:,1:10)))));
shDays = datestr(shUnique)
NDays = length(shUnique)

for i = 1: length(fns)
    i
    load([pnParkerOutFinal2 fns(i).name])
    day(i) = floor(epoch2Matlab(str2num(foo(i,1:10))));
    Hs(i) = runups.RMaxFromBrightest.env.Hs;
    fp(i) = runups.RMaxFromBrightest.env.fp;
    zt(i) = runups.RMaxFromBrightest.env.zt;
    NRunups(i) = length(runups.xMax);
end

uniqueDays = unique(day);
for i = 1: length(uniqueDays)
    ind = find(day == uniqueDays(i));
    HsBar(i) = mean(Hs(ind));
    fpBar(i) = mean(fp(ind));
    N(i) = length(ind);
end

[datestr(uniqueDays') HsBar' fpBar']
[HsBar' fpBar' N']
