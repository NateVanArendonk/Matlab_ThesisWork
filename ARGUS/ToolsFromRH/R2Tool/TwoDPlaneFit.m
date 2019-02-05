function z = TwoDPlaneFit(beta, xy)
%
%   2D planar fit z = -(beta(1)*x + beta(2)*y + beta(3)).  Note z = -h.
%
    z = -(beta(1)*xy(:,1) + beta(2)*xy(:,2) + beta(3));
end