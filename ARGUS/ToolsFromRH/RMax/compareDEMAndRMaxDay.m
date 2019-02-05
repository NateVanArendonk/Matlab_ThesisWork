% compare topo maps to lidar

clear
pnRMax = '/home/ruby/users/holman/research/brightestWork/RMax/OUTPUTRedoDEM/';
fnsR = dir([pnRMax '*RMaxDay.mat']);
pnDEMs = '/home/ruby/users/holman/research/brightestWork/DEMDays/';
fnsD = dir([pnDEMs '*DEMDay.mat']);

% for each overlapping day load both and plot various options
yIndsR = [141:201];
xIndsR = [1: 81];
yIndsD = [1:5:301];
xIndsD = [21:101];
cLims = [-1 3];
dzcLims = [-2 2];
errLims = [-1 1];
for i = 1: length(fnsR)
    load([pnRMax fnsR(i).name])
    xR = RMaxDay.mapped.xm(xIndsR);
    yR = RMaxDay.mapped.ym(yIndsR);
    zR = RMaxDay.mapped.zEst.update.Z(yIndsR, xIndsR);
    zRK = RMaxDay.mapped.zEst.Kalman.Z(yIndsR, xIndsR);
    sigzRK = sqrt(RMaxDay.mapped.zEst.Kalman.P(yIndsR, xIndsR));
    
    load([pnDEMs fnsD(i).name])
    xD = d.x(xIndsD);
    yD = d.y(yIndsD);
    zD = d.DEM(yIndsD,xIndsD);
    dZ = zRK-zD;
    dn(i) = epoch2Matlab(RMaxDay.when);
    
    figure(1); clf; subplot(141)
    subplot(151)
    imagesc(xD, yD, zD)
    axis xy; caxis(cLims); colorbar
    title(['DEM' datestr(dn(i),2)])

    subplot(152)
    imagesc(xR, yR, dZ)
    axis xy; caxis(dzcLims); colorbar
    title('Kalman-DEM')
    
    subplot(153)
    imagesc(xR, yR, zRK)
    axis xy; caxis(cLims); colorbar
    title('Kalman')
    
    subplot(154)
    imagesc(xR,yR, sigzRK)
    axis xy; caxis(errLims); colorbar
    title('Kalman error')
    
    subplot(155)
    imagesc(xR,yR, zR)
    axis xy; caxis(cLims); colorbar
    title('Update')
    
    good = find((~isnan(zR)) & (~isnan(zD)));
    dZ = zR(good)-zD(good);
    
    mndZ(i) = mean(dZ(:));
    stddZ(i) = std(dZ(:));
    rmsdZ(i) = rms(dZ(:));
    range(i) = max(zD(good)) - min(zD(good));
    foo = corrcoef(zR(good), zD(good));
    cor(i) = foo(1,2);
    clear foo
    for j = 1: length(RMaxDay.RMaxes)
        foo(j) = RMaxDay.RMaxes(j).env.Hs;
    end
    Hs(i) = mean(foo);
    bD850(i) = d.beta850;
    bR850(i) = RMaxDay.RMaxes(1).env.beta(171);
    [mndZ(i) stddZ(i) rmsdZ(i) range(i) Hs(i) bD850(i) bR850(i)]
    figure(4); clf
    hist(dZ(~isnan(dZ)))
    figure(5); clf
    lims = [55 100];
    drawProfiles(d,RMaxDay, lims)  
    ylim([0 4])
    %pause
end

% examine some of the bulk error statistics
figure(3); clf
subplot(311);
plot(dn, Hs); datetick('x')
ylabel('H_s (m)');grid on

subplot(312);
plot(dn, mndZ, 'b', dn, stddZ, 'r--', dn, rmsdZ, 'g-.', dn, range, 'k:'); datetick('x')
ylabel('error metrics (m)');grid on; datetick('x')
ylim([-2 5])
a = axis; h=line(a(1:2), [0 0]); set(h, 'color', 'k')
legend('bias', 'std', 'rms', 'range','Location', 'southwest')

subplot(313);
plot(dn, bD850, 'b', dn, bR850, 'r'); datetick('x')
ylabel('slope,y=850');grid on; datetick('x')
legend('DEM', 'R_{max}','Location', 'northwest')
ylim([0 0.15])

% make a figure showing the results for the median rms error as a example
medVal = median(rmsdZ)
[~, medInd] = min(abs(rmsdZ - medVal))
i= medInd;
load([pnRMax fnsR(i).name])
xR = RMaxDay.mapped.xm(xIndsR);
yR = RMaxDay.mapped.ym(yIndsR);
zR = RMaxDay.mapped.zEst.update.Z(yIndsR, xIndsR);
zRK = RMaxDay.mapped.zEst.Kalman.Z(yIndsR, xIndsR);
sigzRK = sqrt(RMaxDay.mapped.zEst.Kalman.P(yIndsR, xIndsR));

load([pnDEMs fnsD(i).name])
xD = d.x(xIndsD);
yD = d.y(yIndsD);
zD = d.DEM(yIndsD,xIndsD);
dZ = zRK-zD;
dn(i) = epoch2Matlab(RMaxDay.when);

figure(1); clf; 
subplot(131)
imagesc(xD, yD, zD)
axis xy; caxis([-1 3]); colorbar
title(['DEM, ' datestr(dn(i),2)])
xlabel('x (m)'); ylabel('y (m)')

subplot(133)
imagesc(xR, yR, dZ)
axis xy; caxis([-1 1]); colorbar
title('Kalman-DEM')
xlabel('x (m)'); ylabel('y (m)')

subplot(132)
imagesc(xR, yR, zRK)
axis xy; caxis([-1 3]); colorbar
title('Kalman')
xlabel('x (m)'); ylabel('y (m)')

