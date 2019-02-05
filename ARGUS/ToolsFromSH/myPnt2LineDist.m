function Dist = myPnt2LineDist(Pnt,Line,WhichAxis)
%__________________________________________________________________________
%DIST = myPnt2LineDist(PNT,LINE,WHICHAXIS)
%
%   Calculate the distance(s) DIST from the point(s) PNT to the line LINE 
%   using an one-step formula.
%
%   PNT must be a Nx2 matrix where N is the number of points over which the
%   distance(s) will be calculated, and the 2 columns represent the X and Y 
%   coordinates of the point(s), on this order.
%
%   LINE must be a vector with two values, where the first specifies the 
%   angular coefficient and the second the intersect of the line.
%
%   WHICHAXIS is an optional argument indicating the axis to which the line
%   coefficients are referred to. Choose 'x' for Xaxis (default) and 'y'
%   for Yaxis (the function is not case-sensitive).
%
%   DIST is a column vector of length N that gives the minimum distances 
%   from each point to the line.
%
%   Rafael Guedes, UoW
%   30/01/2011
%__________________________________________________________________________


%HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH%
%                            Adjust point(s) coordinates so they make sense
if nargin==3
    switch lower(WhichAxis)
        case 'x'
            m = Pnt(:,1);
            n = Pnt(:,2);
        case 'y'
            n = Pnt(:,1);
            m = Pnt(:,2);
        otherwise
            error('WHICHAXIS must be either ''x'' or ''y''')
    end
end


%HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH%
%                           Convert line coefficients to the ax+by+c=0 form
% y = Mx + B   --> Mx - y + B = 0
% ax + by + c = 0

a = Line(1);
b = -1;
c = Line(2);


%HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH%
%                                                 Calculate the distance(s)
Dist = abs(a.*m+b.*n+c)./sqrt(a^2+b^2);





