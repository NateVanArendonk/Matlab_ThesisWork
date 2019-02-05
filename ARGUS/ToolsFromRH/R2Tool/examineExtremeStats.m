% examineExtremeStats
% load all of the edited files and test how RMax compares as well as study
% the relative magnitudes of RMax, R2, R5, etc

clear
pn = '/home/ruby/users/holman/research/RUNUP/discreteRunupTool/OUTPUTParkerFinal2/';
pnRMax = '/home/ruby/users/holman/research/brightestWork/RMax/OUTPUT/';
fns = dir([pn '*mat']);
DBConnect blossom.oce.orst.edu holman '' backup_argus

for i = 1: length(fns)
    load([pn fns(i).name])
    rb(i) = runups.RMaxFromBrightest.x;
    rr(i) = runups.RMax10Min.x;
    dn(i) = epoch2Matlab(runups.when);
    m(i) = runups.RMaxFromBrightest.merit;
    R2(i) = runups.R2.x;
    R5(i) = runups.R5.x;
end
err = rb-rr;

% plot
figure(1); clf
scatter(rr,rb,8,m)
xlabel('x of R_{Max} from runups (m)')
ylabel('x of R_{Max} from brightest (m)')
h = colorbar; set(get(h, 'title'),'string', 'merit')
grid on; caxis([40 140])

figure(2); clf
plot(m,abs(err),'o')    % doesn't appear very useful


% sort by err and show images from worst cases in order
[errSort, ind] = sort(abs(err),2,'descend');

yInd = 171;     % y = 850;
for i = 1:50;
    load([pn fns(ind(i)).name])
    load([runups.RMaxFromBrightest.RMaxDayfn]);
    dns = epoch2Matlab([RMaxDay.RMaxes.when]);
    dt = abs(dns-dn(ind(i)));
    foo = find(dt == min(dt));
    fnb = RMaxDay.RMaxes(foo).fnList(1,:);
    I = imread(FTPPath(fnb));
    figure(4); clf
    imagesc(I); hold on
    g = DBGetImageData(fnb);
    [U,V] = findUV(g.geometry.m, RMaxDay.RMaxes(foo).results(:,1:3));
    plot(U,V,'-')
    
    [Ub,Vb] = findUV(g.geometry.m, [runups.RMax10Min.x runups.ym 0]);
    plot(Ub,Vb, 'r*')
    pause
end

% There are several sources of larger error.  Mostly they correspond to
% confused brightest results.  One measure is how much adjacent estimates
% vary (schizophrenia) or don't work at all (% success).

% exclude these two from stats
keep = 1:length(fns);
err = rb(keep)-rr(keep);
figure(2); hist(abs(err))
N = length(keep);
errSort = sort(abs(err));
err50 = errSort(round(N/2))         % answer was 0.23
err90 = errSort(round(0.9*N))       % answer was 1.44
foo = find(errSort<=1.0, 1, 'last')
foo/length(fns)

% try some rudimentary skill tests
yInds = yInd-10:yInd+10;
for i = 1:length(fns)
    load([pn fns(i).name])
    load([runups.RMaxFromBrightest.RMaxDayfn]);
    dns = epoch2Matlab([RMaxDay.RMaxes.when]);
    dt = abs(dns-dn((i)));
    foo = find(dt == min(dt));
    xbSub = RMaxDay.RMaxes(foo).results(yInds,1);
    goodSub = find(RMaxDay.RMaxes(foo).results(yInds,9));
    fractGood(i) = length(goodSub)/length(yInds);
    mad(i) = mean(abs(diff(xbSub(goodSub))));
    maxDev(i) = max(abs(diff(xbSub(goodSub))));
end

plot(fractGood, abs(err), '*')
plot(mad, abs(err), '*')
plot(mad./fractGood, abs(err), '*')
plot(maxDev, abs(err), '*')

keep = find(fractGood>0.75);
err = rb(keep)-rr(keep);
figure(2); hist(abs(err))
N = length(keep);
errSort = sort(abs(err));
err50 = errSort(round(N/2))         % answer was 0.23
err90 = errSort(round(0.9*N))       % answer was 1.44

