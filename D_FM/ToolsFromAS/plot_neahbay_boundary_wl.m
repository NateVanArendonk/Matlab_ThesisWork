load(['y:\Projects\pscosmos\data\wl\coops_2017\navd88\',...
    'NOAA_9443090_Neah_Bay_navd88.mat'])


    
 addpath('c:\data\tools\CMGmfiles\')
[ntr, jdf]=cmglowpass(wld.WL_VALUE, 360, wld.time,1,'osu');

doy=date2doy(jdf);
minx=floor(min(jdf));
maxx=ceil(max(jdf));

ts=datenum(2017,2,1);
te=datenum(2017,4,1);


figure
lh(1)=plot(wld.time,wld.WL_VALUE);
hold on
lh(2)=plot(jdf,ntr,'linewi',2);

lh(3)=line([minx maxx],[1.059 1.059],'color','k',...
    'linewi',2,'linestyle','--');
set(gca,'xlim',[minx maxx])
yl=ylim;
yp=[yl(2) yl(2) yl(1) yl(1) yl(2)];
xp=[ts te te ts ts];
ph=patch(xp,yp,[0.6 0.6 0.6]);
uistack(ph,'bottom');

datetick('x','mmm','keeplimits')
ylabel('\bf\itWater level (m, NAVD88)')
xlabel('\bf\it2017')

set(gca,'ylim',yl,'layer','top')

leg=legend([lh,ph],'Water level','Non-tidal residual','Station MSL',...
    'Calibration time frame');
