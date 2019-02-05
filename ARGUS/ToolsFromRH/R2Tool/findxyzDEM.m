function xyzDEM = findxyzDEM(xyz0, xyzCam, d)
%   xyzDEM = findxyzDEM(xyz0, xyzCam, dem)
%
%   transform basic RMax xyz0 data (assuming some constant reference level)
%   to true xyz positions based on a supplied DEM with fields DEM.x, DEM.y,
%   DEM.z.

xyzDEM = repmat(nan,size(xyz0));
xyzDEM = pixels2DEM(d.x, d.y, d.z, xyz0, xyzCam);

