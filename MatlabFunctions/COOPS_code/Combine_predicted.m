ln1 = length(one.time);
ln2 = length(two.time);
ln3 = length(three.time);
totl = ln1+ln2+ln3;

% Grab WL
wl1 = one.wl;
wl2 = two.wl;
wl3 = three.wl;

% Grab time
t1 = one.time;
t2 = two.time;
t3 = three.time;

% Create structure with values
pre.wl = [wl1 wl2 wl3];
pre.time = [t1 t2 t3];




