% examineExtremeStats
% load all of the edited files and test how RMax compares as well as study
% the relative magnitudes of RMax, R2, R5, etc.  This file looks only at
% the vertical components of the wave contribution found from the DEM.

clear
pn = '/home/ruby/users/holman/research/RUNUP/discreteRunupTool/OUTPUTParkerFinal3/';
fns = dir([pn '*mat']);
DBConnect blossom.oce.orst.edu holman '' backup_argus

for i = 1: length(fns)
    load([pn fns(i).name])
    dn(i) = epoch2Matlab(runups.when);
    rb(i) = runups.RMaxFromBrightest.x;
    rbz(i) = runups.RMaxFromBrightest.zWave;
    m(i) = runups.RMaxFromBrightest.merit;
    rr(i) = runups.RMax10Min.x;
    rr10z(i) = runups.RMax10Min.zWave;
    rrz(i) = runups.RMax.zWave;
    R2z(i) = runups.R2.zWave;
    R2x(i) = runups.R2.x;
    R5z(i) = runups.R5.zWave;
    R10z(i) = runups.R10.zWave;
    R50z(i) = runups.R50.zWave;
    R2Stock(i) = runups.RMaxFromBrightest.env.R2;
    zt(i) = runups.RMaxFromBrightest.env.zt;
    Hs(i) = runups.RMaxFromBrightest.env.Hs;
    fp(i) = runups.RMaxFromBrightest.env.fp;
    beta(i) = runups.RMaxFromBrightest.env.beta;
    betaLidar(i) = runups.RMaxFromBrightest.env.betaLidar;
    Lo = 9.8/(2*pi*fp(i)*fp(i));
    irr(i) = beta(i)/sqrt(Hs(i)/Lo);
    NMaxes(i) = length(runups.xMax);
end


% plot extreme runup comparison for 10 minute records.
figure(1); clf
scatter(rr,rb,8,m)
xlabel('x of R_{Max} from runups (m)')
ylabel('x of R_{Max} from brightest (m)')
h = colorbar; set(get(h, 'title'),'string', 'merit')
grid on; caxis([40 140])

figure(2); clf
plot(betaLidar, beta, '*')

figure(3); clf
plot(R2Stock, R2z, '*'); grid on; axis([0 3 0 3])
xlabel('R_2 from Stockdon (m)'); ylabel('R_2 from runups (m)')
hold on; 
diss = find(irr<=0.3);
plot(R2Stock(diss), R2z(diss), 'r*')

% R2 versus RMax.
figure(4);clf
plot(R2z,rrz, '*'); axis([0 3.5 0 3.5]); hold on
bad = [36 37 39];
plot(R2z(bad), rrz(bad), 'r*')
xlabel('R_2 from runups (m)'); ylabel('R_{max} from runups (m)')
grid on
% exclude the known bad points
keep = setdiff([1:length(R2z)], bad);
R2z = R2z(keep); rrz = rrz(keep);
good = find((~isnan(R2z)) & (~isnan(rrz)));
m = polyfit(R2z(good),rrz(good),1)  % slope is 1.0056;
k = rrz(good)./R2z(good);
%clf;hist(k, [1:0.02:1.4])
mean(k) %1.071
median(k)    %1.042
[kSort, ind]=sort(k);
kSort(ind(round(0.95*length(k))))   % 95% less than 1.29.

% do other stats 
R5z = R5z(keep);
k = R5z(good)./R2z(good);
mean(k) % 0.908;
median(k) % 0.912
R10z = R10z(keep);
k = R10z(good)./R2z(good);
mean(k) % 0.816
median(k) % 0.830
R50z = R50z(keep);
good = find((~isnan(R50z)) & (~isnan(R2z)));
k = R50z(good)./R2z(good);
mean(k) % 0.566
median(k) % 0.559
rr10z = rr10z(keep);
good = find((~isnan(rr10z)) & (~isnan(R2z)));
k = rr10z(good)./R2z(good);
mean(k) % 1.032
median(k) % 1.002
rr10z = rr10z(keep);
good = find((~isnan(rr10z)) & (~isnan(rrz)));
k = rrz(good)./rr10z(good);
mean(k) % 1.039
median(k) % 1.00


figure(3); clf
plot(




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

