function d = getClosestDEM(e)
%   DEM = getClosestDEM(epoch)
%
%  find the closest LIDAR DEM from a local repository to an epoch time

pn = '/home/ruby/users/holman/research/brightestWork/DEM_20150918_20151031/';
foo = datestr(epoch2Matlab(e),30);
fns = dir([pn foo(1:8) '*']);

if ~isempty(fns)
    for i = 1: length(fns)
        hrs(i) = str2num([fns(i).name(10:11)]);
    end
    [delHr, ind] = min(abs(hrs-str2num(foo(10:11))));
    d.fullPn = [pn, fns(ind).name];
    load(d.fullPn)
    d.x = frameGriddedData.xs;
    d.y = frameGriddedData.as;
    d.z = frameGriddedData.data;
else
    warning(['No DEM for ' datestr(epoch2Matlab(e),2)])
    d.x = [];
    d.y = [];
    d.z = [];
    d.fullPn = [];
end
