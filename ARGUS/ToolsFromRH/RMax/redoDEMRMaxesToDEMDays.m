% read in all RMaxDays.  For every RMaxes, load DEMDay data and fill in 
% DEM corrections.  We should later interp those.  
% Rob Holman, Sept 2018.

clear
addpath('/home/ruby/matlab/CILtest/CIL/argusDB2/');
DBConnect blossom.oce.orst.edu holman '' backup_argus

% find all RMaxDay files
pn = '/home/ruby/users/holman/research/brightestWork/RMax/OUTPUT/';
fns = dir([pn '*RMaxDay.mat']);
pnDEMs = '/home/ruby/users/holman/research/brightestWork/DEMDays/';
fnsD = dir([pnDEMs '*DEMDay.mat']);
inds = 49:92;       % overlapping dates with Lidar
fnsR = fns(inds);

for i = 1: length(fnsR)
    clear RMaxDay d
    load([pn fnsR(i).name])
    DEMpn = [pnDEMs fnsD(i).name];
    load(DEMpn)
    R = RMaxDay.RMaxes;
    DEMpn = [pnDEMs fnsD(i).name];
    for c = 1: length(RMaxDay.camNums)
        camInfo(c).xyz = RMaxDay.xyzCams(c,:);
    end
    for k = 1: length(R)
        R(k).xyzDEM = findxyzDEM(R(k), camInfo, d);
    end
    RMaxDay.RMaxes = R;
    RMaxDay = makeMapsRMaxDay(RMaxDay);      % assemble the full day and save
    eval(['save OUTPUTRedoDEM/' fnsR(i).name ' RMaxDay'])
    %examineRMaxDay(RMaxDay)                 % examine if useful
end

