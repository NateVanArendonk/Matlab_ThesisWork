clearvars
name = 'station_data\PAM ROCKS';
S = load(name);

time = datenum(S.yr, S.mo, S.da, S.hr, S.mn, 30);
clf
plot(time,S.spd)
datetick()

%%
inds = find(S.yr <= 1990);
S.spd(inds) = NaN;
time(inds) = NaN;
clf
plot(time,S.spd)
datetick()


S.alt(inds) = NaN;
S.da(inds) = NaN;
S.gust(inds) = NaN;
S.hr(inds) = NaN;
S.mn(inds) = NaN;
S.mo(inds) = NaN;
S.slp(inds) = NaN;
S.usaf(inds) = NaN;
S.wban(inds) = NaN;
S.wnddir(inds) = NaN;
S.yr(inds) = NaN;
save(name,'-struct','S');


