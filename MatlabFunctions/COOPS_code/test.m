ln1 = length(one.time);
ln2 = length(two.time);
ln3 = length(three.time);
totl = ln1+ln2+ln3;


%tide.wl = NaN(1,totl);



wl1 = one.wl;
wl2 = two.wl;
wl3 = three.wl;

t1 = one.time;
t2 = two.time;
t3 = three.time;


tide.wl = [wl1 wl2 wl3];
tide.time = [t1 t2 t3];

