function xyzEst = findxyzEst(RMaxes, camInfo)
%   RMaxDay = findxyzEst(RMaxDay, camInfo)
%
% Takes the xyz RMax results from findRMax which are based on an assumed z
% value (often 0) and maps them onto the appropriate locations based on an
% estimated z level contained in env.delz (function of y).  This z level
% may come from the Stockdon estimate or from other sources.
% The correction is based on a linear interpolation along the UV rays.

xyz0 = RMaxes.results(:,1:3);       % grad zero level results
c = RMaxes.results(:,4);
cams = unique(c(~isnan(c)));
dz = RMaxes.env.delz;
xyzEst = nan(size(xyz0));
for i = 1:length(cams)              % cycle through cams
    pick = find(c == cams(i));      % pick out RMaxes for this cam
    xyzC = repmat(camInfo(i).xyz, length(pick),1);  % repmat it
    dxyz = (xyz0(pick,:)-xyzC);                 % vector direction
    scale = repmat((dxyz(:,3)+dz(pick))./dxyz(:,3), 1, 3);
    xyzEst(pick,:) = dxyz.*scale + xyzC;   % scale vector to adjust z
end

