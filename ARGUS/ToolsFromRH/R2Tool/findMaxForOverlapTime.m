function r = findMaxForOverlapTime(r, e)
%   runups = findMaxForOverlapTime(runups, et)
%
%  Takes a set of runups data and finds the RMax from the discrete maxima
%  for the period of overlap with a brightest image.  The period of overlap
%  goes from a specified epoch time, et, for 600 seconds.

brightestDur = 600;
ind = find((r.tMax>=e) & (r.tMax<=e+brightestDur));
[RMax10Min, foo] = min(r.xMax(ind));
r.RMax10Min.x = RMax10Min;
r.RMax10Min.t = r.tMax(ind(foo));
r.RMax10Min.ind = ind(foo);
