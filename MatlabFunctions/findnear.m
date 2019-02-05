function [y,pos] = findnear(x,b,ud)
% FINDNEAR			[y,pos] = findnear(x,b,ud)
%
%	FINDNEAR(x,b) finds the element 'y' and it's location
%	'pos' in array 'x' closest to the value 'b'.  If the
%	element is equidistant from the two surrounding
%	elements in array 'x', the user can choose the element
%	before or the element after 'y' using 'ud'.  If
%	'ud' = 1, then the the element before 'y' is choosen;
%	conversely, if 'ud' = 2, then the element after 'y' is
%	choosen.
%
%
%									Curt Storlazzi  03/01


dist = abs(x-b);

%[pos1,pos2] = find(dist == repmat(min(dist(:)),size(dist)));	% The minimum of the distance vector
[ipos] = find(dist == repmat(min(dist(:)),size(dist)));

if length(ipos) >= 2 & ud == 1
	pos = min(min(ipos));
elseif length(ipos) >= 2 & ud == 2
	pos = max(max(ipos));
else
	pos = ipos;
end

%value = A(pos1,pos2);
y = x(pos);