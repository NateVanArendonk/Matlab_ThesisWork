% examine the Sept 18 RMaxDay run for vertical errors

% get list of RMaxDays from OUTPUTRedoetc
clear DEMBar HsBar R2Bar
for i = 1: length(fnsR)
load([pnRMax fnsR(i).name])

yInd = 171;
for j = 1: length(RMaxDay.RMaxes)
    xyz0(j,:) = RMaxDay.RMaxes(j).results(yInd,1:3);
    xyzEst(j,:) = RMaxDay.RMaxes(j).xyzEst(yInd,:);
    Hs(j) = RMaxDay.RMaxes(j).env.Hs;
    beta(j) = RMaxDay.RMaxes(j).env.beta(yInd);
    zt(j) = RMaxDay.RMaxes(j).env.zt;
    R2(j) = RMaxDay.RMaxes(j).env.R2(yInd);
    delz(j) = RMaxDay.RMaxes(j).env.delz(yInd);
    xyzDEM(j,:) = RMaxDay.RMaxes(j).xyzDEM(yInd,:);
end

figure(3); clf
plot(zt)

plot(xyzDEM(:,3) - xyzEst(:,3))

figure(3); clf
plot(xyzDEM(:,1), xyzDEM(:,3))
hold on; plot(xyzEst(:,1), xyzEst(:,3), 'r')

figure(6);clf
plot(zt, 'k'); hold on
plot(xyzDEM(:,3))
plot(xyzEst(:,3),'r')
plot(R2,'g')
drawnow

DEMBar(i) = nanmean(xyzDEM(:,3) - zt')
HsBar(i) = mean(Hs)
R2Bar(i) = nanmean(R2)
pause
end


% mess with Sept 28 as an example.
load('../DEMDays/20150928DEMDay.mat')
figure(7);clf
plot(d.x,d.DEM(151,:))
