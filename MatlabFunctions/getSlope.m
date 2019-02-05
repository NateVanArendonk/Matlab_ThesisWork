function slope = getSlope(X,Y)
%Given X Y pairs, this function calculates slope
%   Rise over run 
%   X should be pair of x's
%   Y should be pair of y's
%   Boom - slope

slope = abs((Y(2) - Y(1))/(X(2) - X(1)));


end

