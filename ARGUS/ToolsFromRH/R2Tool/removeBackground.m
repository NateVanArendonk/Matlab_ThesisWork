function I2 = removeBackground(I, params)
%   I2 = removeBackground(I, params)
%
%  Remove the background intensity profile from a runup time stack using a
%  running filter approach, as used by runupTool

IBack = I;          % easy way to init the array
[N,M] = size(I);
n0 = repmat(params.NFiltDark,1,M);

for i = 2: N    % row-wise calculation
    dI = I(i,:) - IBack(i-1,:);
    n = n0;
    n(dI>0) = params.NFiltLight;
    IBack(i,:) = (n.*IBack(i-1,:) + I(i,:)) ./ (n+1);
end
IBack(1,:) = mean(IBack(2:100,:));
I2 = I-IBack;