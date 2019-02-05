% examine Month of data

clear
pn = '/home/ruby/users/holman/research/brightestWork/RMax/OUTPUT/';
fns = dir([pn '*Oct*RMaxDay*']);

for i = 1: length(fns)
    load([pn fns(i).name])
    examineRMaxDay(RMaxDay)
    pause
end

% look specifically at the Kalman results

for i = 1: length(fns)
    load([pn fns(i).name])
    examineDetailsOfzEstKalman(RMaxDay)
    pause
end

for i = 1: length(fns)
    load([pn fns(i).name])
    showRMaxPicksOnImages(RMaxDay, 2)
end
