dname='y:\Projects\pscosmos\runs\psfm_m6\DFM_OUTPUT_ps_2d\';
fname='ps_2d_his.nc';

info=ncinfo([dname,fname]);
vars=arrayfun(@(x)(x.Name),info.Variables,'un',0)';
stns=cellstr(ncread([dname,fname],'station_name')');

%grab model data for specified stations along transect
branches={'ps_branch1'};
sidx=cellfun(@(x)(strmatch(x,stns)),...
    branches,'un',0);%#ok
istns=sidx{1}(1:44);
gstns=stns(istns);



%grab model data for specified stations
dstns={'NOAA_9443090_Neah_Bay';...
    'NOAA_9444090_Port_Angeles';...
    'NOAA_9444900_Port_Townsend';...
    'NOAA_9447130_Seattle';...
    'NOAA_9446484_Tacoma'};
[~,~,istns2]=intersect(dstns,stns,'stable');
gstns2=stns(istns2);


tstart=datenum(2017,2,1);
tend=datenum(2017,4,1);
% spinup=2; %days

time=ncread([dname,fname],'time');
it=datenum(2017,1,1); %from netcdf attributes
mtime=it+(time./86400);

st=find(mtime>=tstart,1,'first');
ei=find(mtime<=tend,1,'last');
mtime=mtime(st:ei);

lon=cellfun(@(x)(ncread([dname,fname],...
    'station_x_coordinate',...
    [x 1],[1 1])'),num2cell(istns),'un',0);
lat=cellfun(@(x)(ncread([dname,fname],...
    'station_y_coordinate',...
    [x 1],[1 1])'),num2cell(istns),'un',0);

lons=cellfun(@(x)(ncread([dname,fname],...
    'station_x_coordinate',...
    [x 1],[1 1])'),num2cell(istns2),'un',0);
lats=cellfun(@(x)(ncread([dname,fname],...
    'station_y_coordinate',...
    [x 1],[1 1])'),num2cell(istns2),'un',0);

mstruct=gcm(axesm('tranmerc'));
mstruct.falsenorthing=0;
mstruct.falseeasting=500000;
mstruct.geoid=[6378137 0.081819191042815];
mstruct.mapparallels=0;
mstruct.nparallels=1;
mstruct.origin=[0 -123 0];
mstruct.scalefactor=0.9996;
[x,y]=projfwd(mstruct,cell2mat(lat),cell2mat(lon));
[xs,ys]=projfwd(mstruct,cell2mat(lats),cell2mat(lons));
dist=[0;cumsum(sqrt(diff(x).^2+diff(y).^2))];

xc=num2cell(x);
yc=num2cell(y);
xsc=num2cell(xs);
ysc=num2cell(ys);
dc=num2cell(dist);

wl=cellfun(@(x)(ncread([dname,fname],'waterlevel',...
    [x st],[1 numel(mtime)])'),num2cell(istns),'un',0);
wls=cellfun(@(x)(ncread([dname,fname],'waterlevel',...
    [x st],[1 numel(mtime)])'),num2cell(istns2),'un',0);

%repackage
mstr=repmat(struct('station',[],'time',mtime,...
    'x',[],'y',[],'dist',[],'wl',[]),length(istns),1);
[mstr(:).station]=deal(gstns{:});
[mstr(:).lat]=deal(lat{:});
[mstr(:).lon]=deal(lon{:});
[mstr(:).x]=deal(xc{:});
[mstr(:).y]=deal(yc{:});
[mstr(:).dist]=deal(dc{:});
[mstr(:).wl]=deal(wl{:});

mstrs=repmat(struct('station',[],'time',mtime,...
    'x',[],'y',[],'dist',[],'wl',[]),length(istns2),1);
[mstrs(:).station]=deal(gstns2{:});
[mstrs(:).lat]=deal(lats{:});
[mstrs(:).lon]=deal(lons{:});
[mstrs(:).x]=deal(xsc{:});
[mstrs(:).y]=deal(ysc{:});
[mstrs(:).wl]=deal(wls{:});



%grab data from noaa coops website
% cnames=cellfun(@(x)(x(6:12)),gstns2,'un',0);
% opt.start_time=datestr(mtime(1));
% opt.end_time=datestr(mtime(end));
% opt.datum='MSL';
% opt.out_dir='z:\projects\pscosmos\data\';
% wld=cellfun(@(x)(get_coops_erdapp(x,opt)),cnames);
% cnames=cellfun(@(x)(x(6:12)),gstns2,'un',0);
% dval='IOOS_SixMin_Verified_Water_Leve';
% dname='f:\Projects\pscosmos\data\wl\';
% for i=1:length(cnames)
%     cfiles=dir([dname,'*',cnames{i},'*.mat']);
%     fnames={cfiles(:).name}';
%     data=arrayfun(@(x)(x.(dval)),...
%         cellfun(@(x)(load([dname,x])),fnames));
%     %concatenate files if necessary
%     if length(fnames)>1
%         fields=fieldnames(data);
%         dc=struct2cell(data)';
%         dc2=cellfun(@(x)(cell2mat(dc(:,x))),...
%             num2cell(1:length(fields)),'un',0);
%         datar=cell2struct(dc2',fields);
%         
%         %make sure times are unique
%         [~,tidx]=unique(datar.time);
%         wld(i)=structfun(@(x)(x(tidx,:)),datar,'un',0);
%     end
%     
%     %opendap to datenum
%     wld(i).time=(wld(i).time./86400) + 719529;
% end
% wld=wld(:);


dname='y:\Projects\pscosmos\data\wl\coops_2017\navd88\';
fnames=dir([dname,'*.mat']);
corder=[1 2 3 5 4];
wld=arrayfun(@(x)(x.wld),...
    arrayfun(@(x)(load([dname,x.name])),fnames(corder)));
wl_raw=arrayfun(@(x)(double(x.WL_VALUE)),wld,'un',0);
[wld(:).WL_VALUE]=deal(wl_raw{:});

%calculate nearest point 
dstr=ipdm([xs ys],[x y],'subset','nearest',...
    'result','structure');
disto=dist(dstr.columnindex);
% 
% 
% %correct wl data to NAVD88 (same as model bathymetry)
% %some station pages do not include published offsets
% %some found here: https://www.ngs.noaa.gov/Tidal_Elevation/
% %No offset found for Friday Harbor, so removed from analysis for now
% offset=[1.05;...  %noaa datums page (neah)
%     1.1650;...    %noaa datums page (port angeles)
%     1.262;...     %noaa tidal datum page (pid: AI2202, pt townsend)
%     1.3090;...    %noaa datums page (seattle)    
%     1.3560];       %noaa datums page (tacoma)
% 
% 
dlat=arrayfun(@(x)(x.latitude(1)),wld,'un',0);
dlon=arrayfun(@(x)(x.longitude(1)),wld,'un',0);
% 
% wl_msl=arrayfun(@(x)(double(x.WL_VALUE)),wld,'un',0);
% wl_navd=cellfun(@(x,y)(x+y),wl_msl,num2cell(offset),'un',0);
% [wld(:).DATUM]=deal('NAVD88');
% [wld(:).WL_VALUE]=deal(wl_navd{:});

%interp data onto model timeframe
wli=arrayfun(@(x)(interp1(x.time,x.WL_VALUE,mtime)),wld,'un',0);



mm=arrayfun(@(x)(mean(x.wl)),mstr);
ms=arrayfun(@(x)(mean(x.wl)),mstrs);
dm=cellfun(@(x)(nanmean(x)),wli);

vv=arrayfun(@(x)(var(x.wl)),mstr);
vs=arrayfun(@(x)(var(x.wl)),mstrs);
dv=cellfun(@(x)(nanvar(x)),wli);

figure
plot(dist,mm,'ko-')
hold on 
plot(disto,ms,'ro','markerfacecolor','r')
plot(disto,dm,'ks','markerfacecolor','g')





%plot the stations
serverURL = ['https://basemap.nationalmap.gov/arcgis/services/',...
    'USGSTopo/MapServer/WMSServer?request=GetCapabilities&service=WMS'];
info = wmsinfo(serverURL);

latlim=[46.75 48.9];
lonlim=[-125 -122];

imageLength=1024*2;

orthoLayer = info.Layer(1);
 [A, R] = wmsread(orthoLayer, 'Latlim', latlim,...
            'Lonlim', lonlim, ...
            'ImageHeight', imageLength,...
            'ImageWidth', imageLength');

mkrs={'o';'^';'s';'>';'<'};
cols=num2cell(jet(length(mkrs)),2);

latm=cell2mat(lat);
lonm=cell2mat(lon);

fun=@(x)([min(x) max(x)]);
clims=fun(mm);

f=figure;
set(f,'renderer','zbuffer','units','inches',...
    'position',[0 0 9 7],...
    'paperpositionmode','au')

subplot(1,3,1)
axesm mercator
hold on 
geoshow(A,R)
scatterm(cell2mat(lat),...
    cell2mat(lon),10,mm,'filled')


lh=cellfun(@(x,y,z,h)(plotm(x,y,'marker',z,'color','k',...
    'markersize',6,...
    'markerfacecolor',h)),...
    dlat,dlon,mkrs,cols,'un',0);
plotm(latm(dstr.columnindex),lonm(dstr.columnindex),...
    'k+')


setm(gca,'maplonlimit',lonlim,...
    'maplatlimit',latlim,...
    'mlabelloc',1,...
    'mlineloc',1,...
    'plabelloc',1,...
    'plineloc',1,...
    'mlabelparallel','north',...
    'plabelmeridian','west',...
    'parallellab','on',...
    'meridianlab','on',...
    'mlabelround',0,...
    'plabelround',0,...
    'grid','on',...
    'flinewidth',2,...
    'frame','on',...
    'labelformat','none');
colormap(jet)
set(gca,'visible','off',...
    'clim',clims)
lx=lonlim(1)+0.01*diff(lonlim);
ly=latlim(1)+0.05*diff(latlim);

[xr,yr]=projfwd(getm(gca),ly,lx);
rh = scaleruler;
setm(rh,'units','km',...
    'majortick',(0:25:50),...
    'majorticklength',5,...
    'minortick',[0 0],...
    'minorticklength',0,...
    'xloc',xr,...
    'yloc',yr,...
    'color','k',...
    'linewi',2,...
    'fontang','it',...
    'fontweight','b')

ap=get(gca,'position');
ap2=[ap(1)-0.125 ap(2)-0.075 ap(3)+0.19 ap(4)+0.19];
set(gca,'position',ap2);
    
c2=colorbar('horiz');
cp=get(c2,'position');
cp2=cp;
cp2(1)=cp(1)+0.1;
cp2(2)=cp(2)-0.02;
cp2(3)=cp(3)*0.5;
cp2(4)=cp(4)*0.8;
set(c2,'position',cp2);
set(get(c2,'xlabel'),'string',...
    'Mean WL (m, NAVD88)')


ax2(1)=subplot(2,3,2:3);
plot([mstr(:).dist]./1000,mm,'ko-')
hold on 
plot(disto./1000,ms,'ks','markerfacecolor','g')
 plot(disto./1000,dm,'ro','markerfacecolor','r')
 leg=legend('Model centerline','Model stations','Measurements');
set(leg,'location','northwest')
set(gca,'xticklabel',[])
      
    ylabel( 'Mean WL (m, NAVD88)')

ax2(2)= subplot(2,3,5:6);
plot([mstr(:).dist]./1000,vv,'ko-')
hold on 
plot(disto./1000,vs,'ks','markerfacecolor','g')
 plot(disto./1000,dv,'ro','markerfacecolor','r')
ylabel( 'WL variance (m)')
xlabel('\bf\itDistance (km)')

set(ax2,'yaxislocation','right')
 
 
 ap=get(ax2,'position');
 offset=num2cell([-0.075;0.075]);
 ap2=cellfun(@(x,y)([x(1) x(2)+y x(3) x(4)*0.8]),ap,offset,'un',0);
 set(ax2,{'position'},ap2)
 