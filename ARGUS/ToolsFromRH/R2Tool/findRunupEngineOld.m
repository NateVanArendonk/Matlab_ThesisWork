function [xr, rInd] = findRunupEngine(x, I, params)
%   Ir = findRunupEngine(x, IDeBack, params)
%
%  find approx runup from a runup stack whose background has been removed.

[N,M] = size(I);
threshUp = repmat(params.upRushInt,1,M);
% get first estimate, then use this in subsequent
rInd(1) = find(I(1,:)>threshUp, 1, 'first');
for i = 2: N
    thresh = threshUp;
    thresh(rInd(i-1)+1:end) = params.downRushInt;
    foo = find(I(i,:)>thresh, 1, 'first');
    if ~isempty(foo)
        rInd(i) = foo;
    else
        rInd(i) = M;    % if fail, stick at offshore limit
    end
end
xr = x(rInd);