dname='z:\projects\pscosmos\data\bathy\andrew\nc\';
fnames=dir([dname,'*.nc']);

minmax=@(x)([min(x) max(x)]);

lon=arrayfun(@(x)(ncread([dname,x.name],'lon')),fnames,'un',0);
lat=arrayfun(@(x)(ncread([dname,x.name],'lat')),fnames,'un',0);

xl=cellfun(@(x)(minmax(x)),lon,'un',0);
yl=cellfun(@(x)(minmax(x)),lat,'un',0);
px=cellfun(@(x)([x,fliplr(x),x(1)]),xl,'un',0);
py=cellfun(@(x)([x(2) x(2) x(1) x(1) x(2)]),yl,'un',0);

ldb=landboundary('read',['z:\projects\pscosmos\data\bathy\',...
    'salish_sea_shoreline.ldb']);



%plot the stations
serverURL = ['https://basemap.nationalmap.gov/arcgis/services/',...
    'USGSTopo/MapServer/WMSServer?request=GetCapabilities&service=WMS'];
info = wmsinfo(serverURL);

latlim=[46.8 51.8];
lonlim=[-130 -122];

imageLength=1024*2;

orthoLayer = info.Layer(1);
 [A, R] = wmsread(orthoLayer, 'Latlim', latlim,...
            'Lonlim', lonlim, ...
            'ImageHeight', imageLength,...
            'ImageWidth', imageLength');


figure
axesm mercator

cols=get(gca,'colororder');
colc=num2cell(cols(1:length(px),:),2);

geoshow(A,R)
hold on 
sh=cellfun(@(x,y,c)(patchm(y,x,c)),px,py,colc,'un',0);
cellfun(@(x,y)(set(x,'facealpha',0.3,...
    'edgecolor',y,'linewi',2)),sh,colc)
plotm(ldb(:,2),...
    ldb(:,1),'k-')
set(gca,'visible','off')

setm(gca,'maplonlimit',lonlim,...
    'maplatlimit',latlim,...
    'mlabelloc',2,...
    'mlineloc',2,...
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

leg=legend(cat(1,sh{:}),cellfun(@(x)(strtok(x,'.')),...
    {fnames(:).name}','un',0));
set(leg,'location','eastoutside',...
    'box','off','interpreter','none')

inside=@(x,y)(find(x>=y(1) & x<=y(2)));
yind=cellfun(@(x)(inside(x,latlim)),lat,'un',0);
xind=cellfun(@(x)(inside(x,lonlim)),lon,'un',0);

xcount=cellfun(@(x)(numel(x)),xind,'un',0);
ycount=cellfun(@(x)(numel(x)),yind,'un',0);

xs=cellfun(@(x,y)(x(y)),lon,xind,'un',0);
ys=cellfun(@(x,y)(x(y)),lat,yind,'un',0);

files={fnames(:).name}';
zs=cellfun(@(x,y,z,xc,yc)(ncread([dname,x],'Band1',...
            [y(1) z(1)],[xc yc])'),...
            files,xind,yind,...
            xcount,ycount,'un',0);
   
[xm,ym]=cellfun(@(x,y)(meshgrid(x,y)),xs,ys,'un',0);
xyz=cellfun(@(x,y,z)([x(z<=30) y(z<=30) z(z<=30)]),xm,ym,zs,'un',0);

        
           
dout='z:\projects\pscosmos\data\bathy\andrew\nc\xyz\';
cellfun(@(x,n)(saveascii(x,[dout,strtok(n,'.'),'.xyz'],...
   [6 6 2],'\t')),...
    xyz,files)

     
        
%         
% dout='z:\projects\pscosmos\data\bathy\andrew\nc\asc\';
% cellfun(@(n,x,y,z)(arcgridwrite([dout,strtok(n,'.'),'.asc'],...
%     x,y,z','precision','%0.1f','header_precision','0.6f')),...
%     files,xs,ys,zs)
% 
