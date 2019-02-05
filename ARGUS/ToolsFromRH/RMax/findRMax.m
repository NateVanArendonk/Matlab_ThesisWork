function RMax = findRMax(fnList, RMaxDay, camInfo, useBRRatio)
%   RMax = findRMax(fnList, RMaxDay, camInfo, {useBRRatio})
%
% Given a set of synchronous brightest images from each of the available
% cameras, fnList, find the full RMax structure.  RMaxDay contains seed
% info and params.  fnList should contain at least one filename.
% If the optional useBRRatio is true, the analysis is done on the ratio of
% Blue to Red instead of the grayshade image.  This can help isolate bright
% sand from swash.  The default is to use gray (useBRRatio = 0);

if nargin == 3
    useBRRatio = 0;
end
params = RMaxDay.params;
Ny = length(params.ys);
resultsAllCams = nan(Ny, 9);
resultsAllCams(:,9) = 0;        % default to invalid
resultsAllCams(:,2) = params.ys;
resultsAllCams(:,3) = repmat(params.zs,Ny,1);

% set environmental parameters
et = str2num(fnList(1,1:10));
eval(['[zt, Hs, fp] = ' params.envInputProg '(et);']);
RMax.when = et;
RMax.whenDone = matlab2Epoch(now);
RMax.env.zt = zt;
RMax.env.Hs = Hs;
RMax.env.fp = fp;
RMax.env.beta = RMaxDay.beta.Kalman.betaKal;    % use previous day's estimates
eval(['RMax.env.R2 = ' params.findZEstProg '(Hs, fp, RMax.env.beta);']);
RMax.env.delz = RMax.env.R2 + RMax.env.zt;
RMax.env.delx = -RMax.env.delz ./ RMax.env.beta;
RMax.fnList = fnList;
for c = 1: size(fnList,1)        % loop through camera numbers
    fn = fnList(c,:)
    if (fn(1) ~= ' ')            % if blank, there is no image
        % load the gray shade images.  If useBRRatio, create that image
        IRGB = double(imread(FTPPath(fn)));
        if max(max(IRGB(:,:,1))) > 1
            IRGB = IRGB/255;
        end
        if useBRRatio
            R = IRGB(:,:,1);
            B = IRGB(:,:,3);
            I = B./R*255;       % scale to saturated white
            I(isinf(I)) = 10;   % not sure what to put here.
        else
            I = rgb2gray(IRGB);
        end
        results = findBrightEdge(RMaxDay, RMax, I, IRGB, camInfo(c));
        results(:,4) = c;       % temporary way to pass cam nums.
        resultsAllCams = mergeResults(resultsAllCams, results);
        if params.debug.showShoreline % show final shorelines
            figure(c); clf; colormap(gray)
            imagesc(IRGB); hold on; axis off
            good = find((resultsAllCams(:,4) == c) & resultsAllCams(:,9));
            xyz = resultsAllCams(good,1:3);
            [U,V] = findUV(camInfo(c).m, xyz);
            [U,V] = distort(U,V,camInfo(c).cam);
            plot(U,V,'-*')
            xyzSeed(:,2:3) = resultsAllCams(:,2:3); 
            xyzSeed(:,1) = RMaxDay.seed.Kalman.xsKal + RMax.env.delx;
            [U,V] = findUV(camInfo(c).m, xyzSeed);
            [U,V] = distort(U,V, camInfo(c).cam);
            good = find(onScreen(U,V, size(I,2), size(I,1)));
            plot(U(good),V(good),'-*g')
            dn = epoch2Matlab(str2num(RMax.fnList(1,1:10)));
            title(datestr(dn))
            pause(1)
        end
    end
end
RMax.results = resultsAllCams;

