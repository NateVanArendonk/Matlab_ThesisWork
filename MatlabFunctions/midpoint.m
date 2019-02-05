function [xm,ym] = midpoint(x1,y1,x2,y2)
% Calculate midpoint between two points
%   Given two x,y point pairs, calculate the midpoint between them
mid = [(x1+x2)/2, (y1+y2)/2];
xm = mid(1);
ym = mid(2);
end

