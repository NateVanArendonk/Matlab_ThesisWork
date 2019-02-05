base_dir='y:\Projects\pscosmos\runs\roughness_sensitivity\';
mruns=cellfun(@(x)([sprintf('psfm_r%0.0f',x),'\DFM_OUTPUT_ps_2d\']),...
    num2cell(([4 2 1 3 5])'),'un',0);
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


dt=cellfun(@(x,y)(t_tide(x,'start time',mtime{1}(1),...
    'interval',diff(mtime{1}(1:2))*24,...
    'latitude',y)),wli,dlat);


%tidal analysis
const={'M2';'K1'};




%station locations
mamp=cell(length(mruns),length(const));
mpha=cell(length(mruns),length(const));
msamp=cell(length(mruns),length(const));
mspha=cell(length(mruns),length(const));

damp=cell(1,length(const));
dpha=cell(1,length(const));

for j=1:length(const)
            dtidx=arrayfun(@(x)(find(strcmpi(const{j},...
            cellstr(x.name)))),dt,'un',0);
        [dt(:).midx]=deal(dtidx{:});
        damp{j}=arrayfun(@(x)(x.tidecon(x.midx,1)),dt);
        dpha{j}=arrayfun(@(x)(x.tidecon(x.midx,3)),dt);
            bidx=find(abs(diff(dpha{j}))>180)+1;
    dpha{j}(bidx:end)=dpha{j}(bidx:end)+360;
end

for i=1:length(mruns)
    dtidx=arrayfun(@(x)(find(strcmpi(const{j},...
            cellstr(x.name)))),dt,'un',0);
        [dt(:).midx]=deal(dtidx{:});
        damp{j}=arrayfun(@(x)(x.tidecon(x.midx,1)),dt);
        dpha{j}=arrayfun(@(x)(x.tidecon(x.midx,3)),dt);
    
    mstr = cellfun(@(x,z)(t_tide(x,'start time',mtime{1}(1),...
        'interval',diff(mtime{1}(1:2))*24,'latitude',z)),...
        wl{i},lat{i});
    
    mdstr = cellfun(@(x,z)(t_tide(x,'start time',mtime{1}(1),...
        'interval',diff(mtime{1}(1:2))*24,'latitude',z)),...
       wls{i},lats{i});
    
    for j = 1:length(const)
        midx=arrayfun(@(x)(find(strcmpi(const{j},...
            cellstr(x.name)))),mstr,'un',0);
        [mstr(:).midx]=deal(midx{:});
        mamp{i,j}=arrayfun(@(x)(x.tidecon(x.midx,1)),mstr);
        mpha{i,j}=arrayfun(@(x)(x.tidecon(x.midx,3)),mstr);
        
        msidx=arrayfun(@(x)(find(strcmpi(const{j},...
            cellstr(x.name)))),mdstr,'un',0);
        [mdstr(:).midx]=deal(msidx{:});
        msamp{i,j}=arrayfun(@(x)(x.tidecon(x.midx,1)),mdstr);
        mspha{i,j}=arrayfun(@(x)(x.tidecon(x.midx,3)),mdstr);
        
        %fix phases so they don't rotate through 0
        bidx=find(abs(diff(mpha{i,j}))>180)+1;
        mpha{i,j}(bidx:end)=mpha{i,j}(bidx:end)+360;
        
        
        bsidx=find(abs(diff(mspha{i,j}))>180)+1;
        mspha{i,j}(bsidx:end)=mspha{i,j}(bsidx:end)+360;
        
        
    end
end






fname2    = 'ps_2d_fou.nc';
grd         = dflowfm.readNet([base_dir,mruns{3},fname2]);
data  = ncread([base_dir,mruns{3},fname2],'mesh2d_fourier005_amp');
face.mask = true(1,grd.face.FlowElemSize);
tri.mask = face.mask(grd.map3);

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

clims=[0.4 1];

addpath('C:\Program Files\MATLAB\R2017a\toolbox\map\mapproj\')
f=figure;
set(f,'renderer','zbuffer','units','inches',...
    'position',[0 0 9 7],...
    'paperpositionmode','au')

subplot(1,3,1)
axesm mercator
hold on 
geoshow(A,R)
ax1=gca;
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
ap2=[ap(1)-0.085 ap(2)-0.075 ap(3)+0.15 ap(4)+0.15];
set(ax1,'position',ap2);
    
c2=colorbar('horiz');
cp=get(c2,'position');
cp2=cp;
cp2(1)=cp(1)+0.05;
cp2(2)=cp(2)-0.05;
cp2(3)=cp(3)*0.5;
cp2(4)=cp(4)*0.8;
set(c2,'position',cp2);
set(get(c2,'xlabel'),'string',...
    'K1 Amp. (m)')

% legstr={'Uniform salinity = 31';...
%     'Uniform salinity, river input';...
%     'Variable salinity, boundary = 31';...
%     'Variable salinity, boundary = 33';...
%     'Variable salinity, boundary = 35'};
legstr={'Chezy = 55';...
    'Chezy = 60';...
    'Chezy = 62.651';...
    'Chezy = 65';...
    'Chezy = 70'};

ax2(1)=subplot(2,3,2:3);

hold on
lh=cellfun(@(x,y)(plot(dist./1000,x,'o-','color',y,...
    'markerfacecolor',y)),mamp(:,2),cols,'un',0);

leg=legend(cat(1,lh{:}),legstr);
set(leg,'location','northwest',...
    'autoupdate','off',...
    'box','off')

hold on 
cellfun(@(x,y)(plot(disto./1000,x,'+',...
    'color',y)),msamp(:,2),cols)
cellfun(@(x,y,z)(plot(x,y,'marker',z,'color','k',...
    'markersize',8,...
    'markerfacecolor','k')),...
    num2cell(disto./1000),num2cell(damp{:,2}),mkrs,'un',0);

% set(leg,'location','northwest')
set(gca,'box','on','xaxislocation','top')
xlabel('Distance (km)')

      
    ylabel( 'K1 Amp. (m)')
    
    
    ax2(2)= subplot(2,3,5:6);
    hold on
lh=cellfun(@(x,y)(plot(dist./1000,x,'o-','color',y,...
    'markerfacecolor',y)),mpha(:,2),cols,'un',0);


hold on 
cellfun(@(x,y)(plot(disto./1000,x,'+',...
    'color',y)),mspha(:,2),cols)
cellfun(@(x,y,z)(plot(x,y,'marker',z,'color','k',...
    'markersize',8,...
    'markerfacecolor','k')),...
    num2cell(disto./1000),num2cell(dpha{:,2}),mkrs,'un',0);

% set(leg,'location','northwest')

ylabel( 'K1 Pha. (deg)')
set(ax2,'box','on','xaxislocation','top',...
    'yaxislocation','right')

 apc=get(ax2,'position');
 offset=num2cell([-0.075;0.075]);
 apc2=cellfun(@(x,y)([x(1)+0.04 x(2)+y x(3)*0.95 x(4)*0.8]),apc,offset,'un',0);
 set(ax2,{'position'},apc2)
 