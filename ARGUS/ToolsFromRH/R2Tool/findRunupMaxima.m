function [tMax, xMax, maxInd] = findRunupMaxima(t, xr, params)
%   [tMax, xMax, maxInd] = findRunupMaxima(t, xr, params);
%
%  find the runup maxima from a time series of swash runup (x(t)).

dx = diff(xr);
signChange = diff(sign(dx));
if (length(find(abs(signChange)>0)) <= 0)    % no peaks
    tMax = nan;
    xMax = nan;
    maxInd = nan;
else
    % join adjacent +1's into single peak
    plusOnes = find(signChange == 1);
    if plusOnes(end) == length(signChange);
        plusOnes = plusOnes(1:end-1);
    end
    for  i = 1: length(plusOnes)
        if (signChange(plusOnes(i)+1) == 1)
            signChange(plusOnes(i)) = 2;
            signChange(plusOnes(i)+1) = 0;
            i = i+1;
        end
    end

    % first remove short noise spikes
    maxInds = find(signChange == 2);
    minInds = find(signChange == -2);
    minDtInSecs = params.minDt;
    for i = 1: length(maxInds)      % check each for nearby min
        if any(abs(t(maxInds(i))-t(minInds)) < minDtInSecs)
            signChange(maxInds(i)) = 0;
        end
    end
    changeInd = find(signChange ~= 0);
    for i = 1: length(changeInd)-1
        if (signChange(changeInd(i)) == 1) & (signChange(changeInd(i+1)) == -1)
            signChange(changeInd(i)) = -1;      % flip 
        end
        if (signChange(changeInd(i)) == 1) & (sign(dx(changeInd(i+1))) == 1)
            signChange(changeInd(i)) = 0;      % keep only beginning of plateaus 
        end
    end

    maxInd = find(signChange > 0) + 1;
    tMax = t(maxInd); tMax = tMax(:);   % force to column
    xMax = xr(maxInd);
end