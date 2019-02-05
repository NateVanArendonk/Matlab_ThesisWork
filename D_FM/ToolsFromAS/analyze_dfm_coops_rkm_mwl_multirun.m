base_dir='y:\Projects\pscosmos\runs\salinity_sensitivity\';
mruns=cellfun(@(x)([sprintf('psfm_s%0.0f',x),'\DFM_OUTPUT_ps_2d\']),...
    num2cell((1:5)'),'un',0);
fname='ps_2d_his.nc';


info=cellfun(@(x)(ncinfo([base_dir,x,fname])),mruns);
vars=arrayfun(@(y)(arrayfun(@(x)(x.Name),y.Variables,'un',0)'),...
    info,'un',0);
stns=cellfun(@(x)(cellstr(ncread([base_dir,x,fname],'station_name')')),...
    mruns,'un',0);

%grab model data for specified stations along transect
branches={'ps_branch1'};
sidx=cellfun(@(y)(cellfun(@(x)(strmatch(x,y)),...
    branches,'un',0)),stns,'un',0);%#ok
istns=cellfun(@(x)(x{1}(1:44)),sidx,'un',0);
gstns=cellfun(@(x,y)(x(y)),stns,istns,'un',0);



%grab model data for specified stations
dstns={'NOAA_9443090_Neah_Bay';...
    'NOAA_9444090_Port_Angeles';...
    'NOAA_9444900_Port_Townsend';...
    'NOAA_9447130_Seattle';...
    'NOAA_9446484_Tacoma'};
[~,~,istns2]=cellfun(@(x)(intersect(dstns,x,'stable')),stns,'un',0);
gstns2=cellfun(@(x,y)(x(y)),stns,istns2,'un',0);


fname2    = 'ps_2d_fou.nc';
grd         = dflowfm.readNet([base_dir,mruns{3},fname2]);
data  = ncread([base_dir,mruns{3},fname2],'mesh2d_fourier001_mean');
face.mask = true(1,grd.face.FlowElemSize);
tri.mask = face.mask(grd.map3);



time=cellfun(@(x)(ncread([base_dir,x,fname],'time')),mruns,'un',0);

it=datenum(2017,1,1); %from netcdf attributes
mtime=cellfun(@(x)(it+(x./86400)),time,'un',0);

tstart=datenum(2017,2,1);
tend=datenum(2017,4,1);
st=cellfun(@(x)(find(x>=tstart,1,'first')),mtime,'un',0);
ei=cellfun(@(x)(find(x<=tend,1,'last')),mtime,'un',0);
mtime=cellfun(@(x,y,z)(x(y:z)),mtime,st,ei,'un',0);

lon=cellfun(@(z,s,n,i)(cellfun(@(x)(ncread([base_dir,z,fname],...
    'station_x_coordinate',[x 1],[1 1])'),num2cell(i),'un',0)),...
    mruns,st,mtime,istns,'un',0);
lat=cellfun(@(z,s,n,i)(cellfun(@(x)(ncread([base_dir,z,fname],...
    'station_y_coordinate',[x 1],[1 1])'),num2cell(i),'un',0)),...
    mruns,st,mtime,istns,'un',0);


lons=cellfun(@(z,s,n,i)(cellfun(@(x)(ncread([base_dir,z,fname],...
    'station_x_coordinate',[x 1],[1 1])'),num2cell(i),'un',0)),...
    mruns,st,mtime,istns2,'un',0);
lats=cellfun(@(z,s,n,i)(cellfun(@(x)(ncread([base_dir,z,fname],...
    'station_y_coordinate',[x 1],[1 1])'),num2cell(i),'un',0)),...
    mruns,st,mtime,istns2,'un',0);


mstruct=gcm(axesm('tranmerc'));
mstruct.falsenorthing=0;
mstruct.falseeasting=500000;
mstruct.geoid=[6378137 0.081819191042815];
mstruct.mapparallels=0;
mstruct.nparallels=1;
mstruct.origin=[0 -123 0];
mstruct.scalefactor=0.9996;
[x,y]=projfwd(mstruct,cell2mat(lat{1}),cell2mat(lon{1}));
[xs,ys]=projfwd(mstruct,cell2mat(lats{1}),cell2mat(lons{1}));
dist=[0;cumsum(sqrt(diff(x).^2+diff(y).^2))];

xc=num2cell(x);
yc=num2cell(y);
xsc=num2cell(xs);
ysc=num2cell(ys);
dc=num2cell(dist);

wl=cellfun(@(z,s,n,i)(cellfun(@(x)(ncread([base_dir,z,fname],'waterlevel',...
    [x s],[1 numel(n)])'),num2cell(i),'un',0)),mruns,st,mtime,istns,'un',0);
wls=cellfun(@(z,s,n,i)(cellfun(@(x)(ncread([base_dir,z,fname],'waterlevel',...
    [x s],[1 numel(n)])'),num2cell(i),'un',0)),mruns,st,mtime,istns2,'un',0);


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
wli=arrayfun(@(x)(interp1(x.time,x.WL_VALUE,mtime{1})),wld,'un',0);


mwl=cellfun(@(x)(cell2mat(x')),wl,'un',0)';
mm=cellfun(@(x)(mean(x,1)'),mwl,'un',0);
vv=cellfun(@(x)(var(x,1)'),mwl,'un',0);

mwls=cellfun(@(x)(cell2mat(x')),wls,'un',0)';
ms=cellfun(@(x)(mean(x,1)'),mwls,'un',0);
vs=cellfun(@(x)(var(x,1)'),mwls,'un',0);

dm=cellfun(@(x)(nanmean(x)),wli);
dv=cellfun(@(x)(nanvar(x)),wli);



figure
cols=get(gca,'colororder');
cc=num2cell(cols(1:length(mruns),:),2);
hold on
cellfun(@(x,y)(plot(dist,x,'o-','color',y,...
    'markerfacecolor',y)),mm',cc)
hold on 
cellfun(@(x,y)(plot(disto,x,'+',...
    'color',y)),ms',cc)
plot(disto,dm,'kp','markerfacecolor','m',...
    'markersize',10)





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

latm=cell2mat(lat{1});
lonm=cell2mat(lon{1});

fun=@(x)([min(x) max(x)]);
clims=fun(mm{end});

addpath('C:\Program Files\MATLAB\R2017a\toolbox\map\mapproj\')
f=figure;
set(f,'renderer','zbuffer','units','inches',...
    'position',[0 0 9 7],...
    'paperpositionmode','au')

subplot(1,3,1)
axesm mercator
hold on 
geoshow(A,R)
z=0;
[xp,yp]=mfwdtran(grd.node.y,grd.node.x);
p = trisurf(grd.tri(tri.mask,:),xp,yp,zeros(size(xp)),...
    data(grd.map3(tri.mask)));
set(p,'FaceColor','flat','edgecolor','none','facealpha',0.7); 
plotm(cell2mat(lat{1}),...
    cell2mat(lon{1}),'ko','markerfacecolor','k',...
    'markersize',3)


lh=cellfun(@(x,y,z,h)(plotm(x,y,'marker',z,'color','k',...
    'markersize',6,...
    'markerfacecolor','k')),...
    dlat,dlon,mkrs,cols,'un',0);
plotm(latm(dstr.columnindex),lonm(dstr.columnindex),...
    'r+')
[xl,yl]=mfwdtran(latlim,lonlim);


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
    'clim',clims,...
    'xlim',xl,'ylim',yl)
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
ap2=[ap(1)-0.05 ap(2)-0.075 ap(3)+0.19 ap(4)+0.19];
set(gca,'position',ap2);
    
c2=colorbar('horiz');
cp=get(c2,'position');
cp2=cp;
cp2(1)=cp(1)+0.05;
cp2(2)=cp(2)-0.05;
cp2(3)=cp(3)*0.5;
cp2(4)=cp(4)*0.8;
set(c2,'position',cp2);
set(get(c2,'xlabel'),'string',...
    'Mean WL (m, NAVD88)')

legstr={'Uniform salinity = 31';...
    'Uniform salinity, river input';...
    'Variable salinity, boundary = 31';...
    'Variable salinity, boundary = 33';...
    'Variable salinity, boundary = 35'};
ax2(1)=subplot(1,2,2);

hold on
lh=cellfun(@(x,y)(plot(dist./1000,x,'o-','color',y,...
    'markerfacecolor',y)),mm',cols,'un',0);

leg=legend(cat(1,lh{:}),legstr)
set(leg,'location','southoutside',...
    'autoupdate','off',...
    'box','off')

hold on 
cellfun(@(x,y)(plot(disto./1000,x,'+',...
    'color',y)),ms',cols)
cellfun(@(x,y,z)(plot(x,y,'marker',z,'color','k',...
    'markersize',8,...
    'markerfacecolor','k')),...
    num2cell(disto./1000),num2cell(dm),mkrs,'un',0);

% set(leg,'location','northwest')
set(gca,'box','on','xaxislocation','top')
xlabel('Distance (km)')
      
    ylabel( 'Mean WL (m, NAVD88)')
    
    ap=get(gca,'position');
    ap2=ap;
    ap2(1)=ap(1)-0.025;
    ap2(2)=ap(2)-0.010;
    ap2(3)=ap(3)*1.2;
    ap2(4)=ap(4)*0.85;
    set(gca,'position',ap2)
