function r = findRunupsFromStack(r)
%   r = findRunupsFromStack(r);
%
% given RMax from a brightest image, extract values from a runup stack.

fn = findArgusImages(r.when,'argus02b','x',['runup' num2str(r.ym)], 'mat');
r = analyzeRunups(r, fn);
