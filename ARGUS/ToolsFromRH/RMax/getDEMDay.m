function d = getDEMDay(e)
%   DEM = getDEMDay(epoch)
%
%  find the closest DEMDay from a local repository to an epoch time

pn = '/home/ruby/users/holman/research/brightestWork/DEMDays/';
fns = dir([pn '*mat']);
foo = datestr(epoch2Matlab(e),26);
nm = foo([1:4 6:7 9:10]);
load([pn nm 'DEMDay.mat'])

if isempty(d)
    warning(['No DEM for ' datestr(epoch2Matlab(e),2)])
    d.x = [];
    d.y = [];
    d.z = [];
    d.fullPn = [];
end
