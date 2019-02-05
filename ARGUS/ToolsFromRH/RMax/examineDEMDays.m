%   createDEMDay
%  read in each day of DEM data and average them to get a robust product.

clear
pn = '/home/ruby/users/holman/research/brightestWork/DEMDays/';
fns = dir([pn '*mat']);
load([pn fns(1).name]);
d1 = d.DEM;
for i = 2: length(fns)
    load([pn fns(i).name])
    d2 = d.DEM;
    figure(1); clf
    imagesc(d.x, d.y, d1);
    figure(2); clf
    imagesc(d.x,d.y, d2);
    title(fns(i).name)
    figure(3); clf
    imagesc(d.x,d.y, d2-d1)
    caxis([-0.5 0.5]); colorbar
    pause
    d1=d2;
end


