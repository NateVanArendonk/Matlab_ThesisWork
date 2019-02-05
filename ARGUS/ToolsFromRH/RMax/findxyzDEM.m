function xyzDEM = findxyzDEM(RMaxes, camInfo, d)
%   RMaxes = findxyzDEM(RMaxes, camInfo, DEM)
%
%   transform basic RMax xyz0 data (assuming some constant reference level)
%   to true xyz positions based on a supplied DEM with fields DEM.x, DEM.y,
%   DEM.z.

xyz0 = RMaxes.results(:,1:3);
xyzDEM = repmat(nan,size(xyz0));
c = RMaxes.results(:,4);
camList = unique(c(find(~isnan(c))));
% loop through the cameras.
for i = 1:length(camList)
    ind = find(c==camList(i));
    xyzc = camInfo(i).xyz;
    [xyzDEM(ind,:)] = pixels2DEM(d.x, d.y, d.DEM, xyz0(ind,:), xyzc);
end
