function resultsAllCams = mergeResults(resultsAllCams, results)
%   resultsAllCams = mergeResults(resultsAllCams, results);
%
% Merge a new camera-specific set of shorelines, results, with the previous
% cumulative results, resultAllCams.  This is done by finding valid values
% in results (based on x being not a nan), then comparing the perceived
% value of those estimates with those already available in resultsAllCams.
% The value here is pretty ad hoc and likely will change with experience.

good = find(results(:,9));      % valid new estimates
preferNew = results(good,5) .* results(good,6); % R2 times amplitude
preferOld = resultsAllCams(good,5) .* resultsAllCams(good,6);
preferOld(isnan(preferOld)) = 0;
keep = preferNew > preferOld;
resultsAllCams(good(keep),:) = results(good(keep),:);

