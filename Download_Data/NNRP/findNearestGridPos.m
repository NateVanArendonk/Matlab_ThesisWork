function [ i, j ] = findNearestGridPos( grid_lat, grid_lon, pos_lat, pos_lon, min_dist )
%[ i, j ] = findNearestGridPos( grid_lat, grid_lon, grid_i, grid_j, pos_lat, pos_lon )
%   Vectors of grid lat/lon and indices i,j are used
%   position lat/lon are used to find nearest grid locs and returne in i,j


%Constant
radius_earth = 6.371e6; %[m]

% %vectorize
% grid_lat = grid_lat(:);
% grid_lon = grid_lon(:);

% Pythag dist in meters (great circle)
pos_lon_rep = pos_lon*ones(size(grid_lon));
pos_lat_rep = pos_lat*ones(size(grid_lon));
pt_dist = distance('gc',pos_lat_rep, pos_lon_rep, grid_lat, grid_lon)*pi/180*radius_earth;

[mymin, I] = min(pt_dist(:));
[i,j] = ind2sub(size(pt_dist),I);

% Check if too far away
if mymin > min_dist
    error('Lat,Lon position of extraction point is too far from nearest grid cell')
end


end





