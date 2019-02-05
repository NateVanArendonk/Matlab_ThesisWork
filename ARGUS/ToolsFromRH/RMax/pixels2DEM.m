function [xyzOut, dz] = pixels2DEM(x, y, DEM, xyz, xyzCam)
%   [xyzOut, dz] = pixels2DEM(x, y, DEM, xyz, xyzCam)
%
%   given a Digital Elevation Model, DEM(x,y), map the locations of a set
%   of pixels whose 'guess' location was xyz onto the DEM.  xyzCam is the
%   known camera location, xyzOut is the result.  dz is the error of the
%   iterative search after NInterations iterations and can be checked if you are
%   unsure but can otherwise be ignored.

NIterations = 6;              % iterations in search
N = size(xyz,1);
zCam = repmat(xyzCam(3),N,1);
for i = 1:NIterations   
    xyzC = repmat(xyzCam,N,1);
    r = xyz-xyzC;
    dz = interp2(x, y, DEM, xyz(:,1), xyz(:,2)) - xyz(:,3);
    s = (zCam-xyz(:,3)-dz)./(zCam-xyz(:,3));
    xyz = xyzC + r.*repmat(s,1,3);
end
xyzOut = xyz;