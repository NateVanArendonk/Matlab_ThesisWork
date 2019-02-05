%   createDEMDay
%  read in each day of DEM data and average them to get a robust product.

clear
pn = '/home/ruby/users/holman/research/brightestWork/DEM_20150918_20151031/';
pnOut = '/home/ruby/users/holman/research/brightestWork/DEMDays/';
fns = dir([pn '*mat']);
fnStr = [fns.name];
names = reshape(fnStr',40,length(fnStr)/40)';
mn = str2num(names(:,5:6));
day = str2num(names(:,7:8));
uniqueMn = unique(mn);
for i = 1: length(uniqueMn)
    fnsp = fns(mn==uniqueMn(i));
    dayp = day(mn==uniqueMn(i));
    uniqueDay = unique(dayp);
    for j = 1: length(uniqueDay)
        pick = find(dayp == uniqueDay(j));
        for k = 1: length(pick)
            load([pn fnsp(pick(k)).name])
            x = frameGriddedData.xs;
            y = frameGriddedData.as;
            DEMs(:,:,k) = frameGriddedData.data;
        end
        DEM = nanmean(DEMs,3);
        figure(1); clf
        imagesc(x,y,DEM)
        figure(2); clf
        sig = nanstd(DEMs,1,3);
        imagesc(x,y,sig)
        colorbar
        title(fnsp(pick(k)).name)
        figure(3); clf
        plot(x,sig(151,:)); ylabel('std (m)'); title('y=850')
        ylim([0 0.3])
        % finally do a planar fit at y = 850, index = 151.  Do +/-100
        [X,Y] = meshgrid(x,y);
        good = find((X>65) &  (X<75) & (Y>=850-100) & (Y<=850+100));
        xs = X(good); ys = Y(good); zs = DEM(good);
        good = find(~isnan(zs));
        zs = zs(good); xs = xs(good); ys = ys(good)-850;
        beta0 = [0.1 0 min(xs)];
        [betaAll, Resid, J, covB] = nlinfit([xs(:) ys(:)], zs, @TwoDPlaneFit, beta0);
        d.x = x; d.y = y; d.DEM = DEM; d.beta850 = betaAll(1)
        outName = [fnsp(pick(k)).name(1:8) 'DEMDay'];
        figure(3); clf
        plot(d.x,d.DEM(151,:)); axis([40 80 0 8]); 
        x = 65:80; y = repmat(0, size(x));
        zPred = TwoDPlaneFit(betaAll, [x' y']);
        hold on; plot(x,zPred, 'k')
        %pause
        eval(['save ' pnOut outName ' d'])
    end
end


