function [ i, j ] = findNearestGridPos( grid_lat, grid_lon, grid_i, grid_j, pos_lat, pos_lon, min_dist )
%[ i, j ] = findNearestGridPos( grid_lat, grid_lon, grid_i, grid_j, pos_lat, pos_lon )
%   Vectors of grid lat/lon and indices i,j are used
%   position lat/lon are used to find nearest grid locs and returne in i,j
%   Nearest 4 locations are identified, reported in order of increasing
%   distance

%Constant
radius_earth = 6.371e6;

% Pythag dist in meters (great circle)
pos_lon_rep = pos_lon*ones(size(grid_lon));
pos_lat_rep = pos_lat*ones(size(grid_lon));
pt_dist = distance('gc',pos_lat_rep, pos_lon_rep, grid_lat, grid_lon)*pi/180*radius_earth;

% Find 4 nearest grid cells, closest first
for nn = 1:4    
    % minimum and index
    [mymin, I] = min(pt_dist);
    
    % Check if too far away
    if mymin > min_dist
        error('Lat,Lon position of extraction point is too far from nearest grid cell')
    end
    
    % Extract i,j values
    i(nn) = grid_i(I);
    j(nn) = grid_j(I);
    
    % Remove index and start again
    pt_dist(I) = [];
    grid_i(I) = [];
    grid_j(I) = [];  
end


end

