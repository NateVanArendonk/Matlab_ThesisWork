function drawProfiles(d,RMaxDay, lims)

DEMInd = 151;
RInd = 171;

xD = d.x; zD = d.DEM(DEMInd,:);
bD = d.beta850;
plot(xD, zD)
hold on; xlim(lims); grid on

xR = RMaxDay.mapped.xm;
zR = RMaxDay.mapped.zEst.Kalman.Z(RInd,:);
bR = RMaxDay.RMaxes(1).env.beta(RInd);
plot(xR, zR, 'r')

zRu = RMaxDay.mapped.zEst.update.Z(RInd,:);
plot(xR, zRu, 'g')
legend('DEM', 'zEstKal', 'zEstUp')

% draw plane slopes for DEM and RMax just for sanity check
good = find(~isnan(zD));
zDg = zD(good); xDg = xD(good);
foo = find(zDg <= median(zDg), 1, 'first');
x0 = xDg(foo); z0 = zDg(foo);
z1 = z0-bD*(lims(1)-x0);
z2 = z0-bD*(lims(2)-x0);
line(lims, [z1 z2])

good = find(~isnan(zR));
if any(good)
    zRg = zR(good); xRg = xR(good);
    foo = find(zRg <= median(zRg), 1, 'first');
    x0 = xRg(foo); z0 = zRg(foo);
    z1 = z0-bR*(lims(1)-x0);
    z2 = z0-bR*(lims(2)-x0);
    h=line(lims, [z1 z2]);
    set(h, 'color', 'r')
end
