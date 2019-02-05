% show runups on brightest.
% all from info in runups structure
clear
pn = '/home/ruby/users/holman/research/RUNUP/discreteRunupTool/OUTPUTParkerFinal3/';
fns = dir([pn '*mat']);
DBConnect blossom.oce.orst.edu holman '' backup_argus

for i = 1: length(fns)
    load([pn fns(i).name])
    fn = findArgusImages(runups.when, 'argus02b', 1, 'bright', 'jpg');
    I = imread(FTPPath(fn));
    figure(4); clf
    imagesc(I); hold on; title(safeString(fn))
    g = DBGetImageData(fn);
    [U,V] = findUV(g.geometry.m, runups.RMax.xyzDEM);
    plot(U,V,'r*')
    [U,V] = findUV(g.geometry.m, runups.RMaxFromBrightest.xyzDEM);
    plot(U,V,'ro')
    [U,V] = findUV(g.geometry.m, [runups.RMax.x runups.ym 0]);
    plot(U,V,'go')

    % load the DEM
    DEM = getLowTideDEM(runups.when);
    figure(2); clf
    imagesc(DEM.x, DEM.y, DEM.z); colorbar; axis xy
    xyz = runups.RMax.xyzDEM;
    hold on; plot(xyz(1), xyz(2), 'r*')
    dy = DEM.y - xyz(2);
    [~,ind] = min(abs(dy));
    figure(3); clf
    plot(DEM.x,DEM.z(ind,:)); hold on
    plot(xyz(1), xyz(3), 'r*');  grid on
    pause
end

