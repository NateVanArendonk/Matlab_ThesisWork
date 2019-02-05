function showRMaxPicksOnImages(RMaxDay, dec)
%   showRMaxPicksOnImages(RMaxDay, {decimate})
%
%  shows images from all cameras on figures with corresponding numbers for
%  every collection in a day RMaxDay.  This provides a visual assessment of
%  how well the algorithm is doing.  decimate allows you to skip frames, so
%  decimate=2 would show every other frame, dec=3 would show every third
%  frame, etc.

if nargin ==1
    dec = 1;
end
% grab cam info for first image of the day, assume it is unchanging
n = parseFilename(RMaxDay.RMaxes(1).fnList(1,:));
camInfo = findCamInfo(n.station, str2num(n.time));
for i = 1: dec: length(RMaxDay.RMaxes)
    r = RMaxDay.RMaxes(i).results;
    fnList = RMaxDay.RMaxes(i).fnList;
    for j = 1: size(fnList,1)
        figure(j); clf
        if (fnList(j,1) ~= ' ')
            I = imread(FTPPath(fnList(j,:)));
            imagesc(I); hold on
            [U,V] = findUV(camInfo(j).m, [r(:,1:3)]);
            [U,V] = distort(U,V,camInfo(j).cam);
            good = find(onScreen(U,V, size(I,2), size(I,1)));
            plot(U(good),V(good),'*-')
            title(datestr(epoch2Matlab(RMaxDay.RMaxes(i).when)))
            xList = RMaxDay.seed.Kalman.xsKal+RMaxDay.RMaxes(i).env.delx;
            [U,V] = findUV(camInfo(j).m, [xList r(:,2:3)]);
            [U,V] = distort(U,V,camInfo(j).cam);
            good = find(onScreen(U,V, size(I,2), size(I,1)));
            plot(U(good),V(good), 'g+-')
        end
    end
    pause
end
