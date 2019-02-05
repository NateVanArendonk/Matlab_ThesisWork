function [xr, rInd, IMax] = findRunupEngine(x, I, params)
%   [xr, rInd, IMax] = findRunupEngine(x, IDeBack, params)
%
%  find runup from a runup stack whose background has been removed.
%  Returns rInd (stack index of runup), xr (associated x location), and
%  IMax, the max intensity above threshold, averaged over all samples

[N,M] = size(I);
IMax = mean(max(I'));       % find typical runup edge brightness
if IMax < params.minIMax    % not enough signal, abort
    rInd = repmat(M, N, 1);
else                        % adequate signal
    threshUp = repmat(params.upRushFraction*IMax,1,M);
    threshDown = params.downRushFraction*IMax;
    % get first estimate, then use this in subsequent
    try
        rInd(1) = find(I(1,:)>threshUp, 1, 'first');
        for i = 2: N
            thresh = threshUp;
            thresh(rInd(i-1)+1:end) = threshDown;
            foo = find(I(i,:)>thresh, 1, 'first');
            if ~isempty(foo)
                rInd(i) = foo;
            else
                rInd(i) = M;    % if fail, stick at offshore limit
            end
        end
    catch
        rInd=repmat(M,N,1);
    end
end
rInd = rInd(:);
xr = x(rInd);