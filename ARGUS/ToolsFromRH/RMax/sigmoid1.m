function y = sigmoid1(beta, x)
%    y = sigmoid1(beta, x)
%
% compute the 1D sigmoid function for independent values x and parameters,
% beta, a 4 by 1 vector.  Beta = [baseLevel amp rate xCenter] such that
%   y = baseLevel + amp / (1 + exp(rate*(x-xCenter)));

y = beta(1) + beta(2)./(1+exp(-beta(3)*(x-beta(4))));
