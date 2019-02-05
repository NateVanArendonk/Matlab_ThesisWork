function results = findBrightEdge(RMaxDay, RMax, I, IRGB, camInfo)
%       results = findBrightEdge(RMaxDay, I, IRGB, camInfo);
%
% extract a set of maximum runup locations from a brightest image, I, with
% geometry, m.  xyz0 is the initial guess at the shoreline, usually created
% from a previous image and is an N by 3 matrix.  The new analysis uses
% that set of y-locations and z values for the search while the x value is 
% used as a seed for the search.  IRGB is the original images, only passed
% for debugging purposes for transect plots.
% The output is returned as an N by 9 matrix, results, with columns
% corresponding to:
%   xRMax, y, z0, camera, R-squared, amplitude, lambda, I0 and QC, with the 
% first and cols 5-8 described by the sigmoid1 function.  The ninth column
% is a binary assessment of valid or not, judged by the function
% brightestQC.
% The camera column is filled in by the calling routine.
%
% The analysis is done in three steps.  First, data are extracted along 
% transects, then a median filter is applied to remove local spikes.  
% Then,  the strongest edge is found using a long Sobel
% filter.  In step 3, the rising end in this sub-region is fit to a sigmoid
% function, returning the above results.
%
% Specific parameters for the analysis are specified by params. 
% Params must have fields 
%   xRes, the desired cross-shore spatial resolution for transit extraction
%   halfWin, search window half width around initial guess location
%   halfSobelWIn, half width of sobel filter for step 1.
%   halfSubWin, reduced window for step 2 sigmoid1 fit
%   guessLambda, the seed guess for lambda in sigmoid1 call

params = RMaxDay.params;
y = params.ys;
Ny = length(y);
z = params.zs;
if length(z) == 1   % expand scalar
    z = repmat(z, Ny, 1);
end
[NV, NU] = size(I);

NS = round(params.halfSobelWidth/params.xRes);    % Sobel filt in pixels
NSub = round(params.halfSubWin/params.xRes);      % nlinfit window size (pixels)
beta = nan(Ny,4);       % initialize nlinfit params
results = nan(Ny,9);    % initialize total results

% instantaneous shoreline is landward for ref shoreline
xt = RMaxDay.seed.Kalman.xsKal + RMax.env.delx;
xyz0 = [xt y z];
results(:,2:3) = xyz0(:,2:3);   % fill in y, z.

for i = 1: Ny
    x = [xt(i)-params.halfWin: params.xRes: xt(i)+params.halfWin]';
    Nx = length(x);
    xyz = [x repmat(xyz0(i,2:3), Nx, 1)];
    [U,V] = findUV(camInfo.m, xyz);
    [U,V] = distort(U,V,camInfo.cam);
    good = find(onScreen(U, V, NU, NV));
    if length(good) > params.minGood*Nx
        ITrans0 = interp2(1:NU, 1:NV, I, U(good), V(good));
        x0 = x(good);        % potentially reduce to onScreen only
        N = round(params.medFiltLength/params.xRes);
        ITrans = medfilt1(ITrans0, N);
        keep = [ceil(N/2): length(x0)-ceil(N/2)-1];
        x = x0(keep);
        ITrans = ITrans(keep);
        Nx = length(x);
        filt = [ones(NS, 1); 0; -ones(NS, 1)]; % for advancing bright
        IFilt = conv(ITrans, filt, 'valid');
        [~,ind] = max(IFilt);
        ind = ind + (length(filt)-1)/2 + 1;
        % now narrow the search to fine tune and to return CIs
        indSub = [max(ind-NSub,1): min(ind+NSub,Nx)];
        xSub = x(indSub);
        ISub = ITrans(indSub);
        beta0 = [min(ISub) max(ISub)-min(ISub) 1/params.guessLambda x(ind)];
        warning('off')
        beta = nlinfit(xSub, ISub, 'sigmoid1', beta0);
        warning('on')
        results(i,[1 6 7 8]) = beta([4 2 3 1]);
        IPred = sigmoid1(beta, xSub);
        dI = IPred-ISub;
        results(i,5) = 1 - sum(dI.^2)/sum((ISub-mean(ISub)).^2);
        if params.debug.showTransects
            figure(10); clf; colormap(gray)
            imagesc(IRGB); hold on
            plot(U,V);              % transect
            dn = epoch2Matlab(str2num(RMax.fnList(1,1:10)));
            [U,V] = findUV(camInfo.m, results(i,1:3));
            [U,V] = distort(U,V,camInfo.cam);
            h=plot(U,V, '*');
            set(h, 'markersize',10)
            axis off; title(['y = 850 m, ' datestr(dn)])
            [U,V] = findUV(camInfo.m, [xt(i) results(i,2:3)]);
            [U,V] = distort(U,V,camInfo.cam);
            h=plot(U,V, 'g*');
            set(h, 'markersize',10)
            figure(11); clf
            plot(x0,ITrans0); hold on
            plot(x,ITrans, 'k');
            plot(xSub, IPred, 'r')
            xlabel('x (m)'); ylabel('I (arb)');
            title(['y = ' num2str(xyz0(i,2))])
            hold on; plot(results(i,1), sigmoid1(beta, results(i,1)), 'r*')
            results(i,5:9)
            legend('I-orig', 'I-medFilt', 'I-Fit', 'location', 'southeast')
            pause
        end
    end
end
results(isnan(results(:,5)),5) = 0;       % force R2 to zero for nan cases.
results(:,9) = brightestQC(results, params.QCLims);

