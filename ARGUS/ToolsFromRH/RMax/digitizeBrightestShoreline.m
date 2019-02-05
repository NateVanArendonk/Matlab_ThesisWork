function xyz0 = digitizeBrightestShoreline(fnList, delx, params, camInfo)
%   xyz0 = digitizeBrightestShoreline(fnList, delx, params, camInfo);
%
% given a list of brightest filenames for a particular instant in time,
% digitize the shoreline contour, adjusted for anticipated corrections to
% z=0 (correction is delx).  params contains useful parameters.

xyz = [];
for i = 1: size(fnList,1)        % loop through cameras
    fn = fnList(i,:);
    I = double(rgb2gray(imread(FTPPath(fn))));
    figure(1); clf; colormap(gray)
    imagesc(I)
    disp('Click shoreline points in Fig 1.  CR to exit')
    UV = ginput;
    UVu = undistort(UV, camInfo(i).cam);
    xyzp = findXYZ(camInfo(i).m, UVu, params.zs, 3);
    xyz = [xyz; xyzp];
end
xys = sortrows(xyz,2);
x0 = interp1(xys(:,2), xys(:,1), params.ys);
x0 = x0 + delx;             % correct to z = 0 shoreline
xyz0 = [x0, params.ys, repmat(params.zs, size(x0))];
