% compare topo maps to lidar

clear
pnRMax = '/home/ruby/users/holman/research/brightestWork/RMax/OUTPUTRedoDEM/';
fnsR = dir([pnRMax '*RMaxDay.mat']);
pnDEMs = '/home/ruby/users/holman/research/brightestWork/DEMDays/';
fnsD = dir([pnDEMs '*DEMDay.mat']);

% for each overlapping day load both and plot various options
yInd = 171;         % pick out y = 850;
k = 1;
for i = 1: length(fnsR)  
    i
    load([pnRMax fnsR(i).name])
    for j = 1: length(RMaxDay.RMaxes)
        x0(k) = RMaxDay.RMaxes(j).results(yInd,1);
        xE(k) = RMaxDay.RMaxes(j).xyzEst(yInd,1);
        xD(k) = RMaxDay.RMaxes(j).xyzDEM(yInd,1);
        t(k) = RMaxDay.RMaxes(j).when;
        zt(k) = RMaxDay.RMaxes(j).env.zt;
        k = k+1;
    end
end

figure(1); clf
plot(x0, xD, '.')
xlabel('x0 (m)'); ylabel('xDEM (m)');
grid on

figure(2); clf
dx = x0-xD;
hist(x0-xD)
nanmean(dx)
max(dx)

max(xD)
min(xD)
max(xD)-min(xD)



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
    
    figure(1); clf; subplot(141)
    subplot(151)
    imagesc(xD, yD, zD)
    axis xy; caxis(cLims); colorbar
    title(['DEM' datestr(epoch2Matlab(RMaxDay.when),2)])

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
    [mndZ(i) stddZ(i) range(i) Hs(i) bD850(i) bR850(i)]
    figure(4); clf
    hist(dZ(~isnan(dZ)))
    figure(5); clf
    lims = [55 100];
    drawProfiles(d,RMaxDay, lims)  
    ylim([0 4])
    pause
end
