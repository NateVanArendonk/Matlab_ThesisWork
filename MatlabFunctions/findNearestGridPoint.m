function [ x_ind, y_ind ] = findNearestGridPoint( x_grid, y_grid, x_point, y_point )
%[ x_ind, y_ind ] = findNearestGridPoint( x_grid, y_grid, x_point, y_point )

dist = (x_grid-x_point).^2 + (y_grid-y_point).^2;

my_min = min(dist(:));

[x_ind, y_ind] = find(my_min==dist);


end