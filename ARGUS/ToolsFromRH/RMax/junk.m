% go through results to estimate how often we get flakey results

clear
addpath('/home/ruby/matlab/CILtest/CIL/argusDB2/');
DBConnect blossom.oce.orst.edu holman '' backup_argus

pn = '/home/ruby/users/holman/research/brightestWork/RMax/OUTPUTRedoDEM/';
fns = dir([pn '*RMaxDay*mat']);

for i = 1: length(fns)
    load([pn fns(i).name])
    showRMaxPickOnImagesSingleCamYPick(RMaxDay,1,1,850);
end
