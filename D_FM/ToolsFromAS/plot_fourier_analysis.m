dname='y:\Projects\pscosmos\runs\salinity_sensitivity\psfm_s5\DFM_OUTPUT_ps_2d\';
fname    = 'ps_2d_fou.nc';
grd         = dflowfm.readNet([dname,fname]);

vars={'mesh2d_fourier001_mean';...
    'mesh2d_fourier002_amp';...
    'mesh2d_fourier002_phs';...
    'mesh2d_fourier005_amp';...
    'mesh2d_fourier005_phs'};


data  = cellfun(@(x)(ncread([dname,fname],x)),vars,'un',0);
face.mask = true(1,grd.face.FlowElemSize);
tri.mask = face.mask(grd.map3);

ldb=landboundary('read',['y:\Projects\pscosmos\runs\psfm_m1\',...
    'salish_sea_shoreline.ldb']);

bname='y:\Projects\pscosmos\data\bathy\andrew\nc\gebco\';
bfile='GEBCO_2014_2D_-130.0_46.8_-122.0_51.8.nc';
lat=ncread([bname,bfile],'lat');
lon=ncread([bname,bfile],'lon');
z=ncread([bname,bfile],'elevation');
xd=vdist(lat(1),lon(1),lat(1),lon(end));
yd=vdist(lat(1),lon(1),lat(end),lon(end));
xi=linspace(0,xd,numel(lon));
yi=linspace(0,yd,numel(lat));
hs=real2rgb(hillshade(z',xi,yi,'zfactor',2),'gray',[0 255]);


fun=@(x)([min(x) max(x)]);
lonlim=[-129.22892620896	-122.084659449228];
latlim=[46.9675525931304	51.6220608123692];


[lonm,latm]=meshgrid(lon,lat);

pvals={[1 1.6],'Mean WL (m, NAVD88)';...
    [0.2 1.6],'M2 Amplitude (m)';...
    [0.4 0.8],'K1 Amplitude (m)'};
aparam=[1 2 4];
pidx=1;    


addpath('C:\Program Files\MATLAB\R2017a\toolbox\map\mapproj\')
f=figure;
set(f,'renderer','zbuffer','units','inches',...
    'position',[0 0 9 7],...
    'paperpositionmode','au')

axesm('mercator')
hold on 
[x,y]=mfwdtran(latm,lonm);
[xp,yp]=mfwdtran(grd.node.y,grd.node.x);
[xl,yl]=mfwdtran(latlim,lonlim);

sh=surf(x,y,zeros(size(x)),hs);
set(sh,'edgecolor','none')

hold on

z=0;
p = trisurf(grd.tri(tri.mask,:),xp,yp,zeros(size(xp)),...
    data{aparam(pidx)}(grd.map3(tri.mask)));
set(p,'FaceColor','flat','edgecolor','none','facealpha',0.7); 

plotm(ldb(:,2),ldb(:,1),'k-')


setm(gca,'maplonlimit',lonlim,...
    'maplatlimit',latlim,...
    'mlabelloc',2,...
    'mlineloc',2,...
    'plabelloc',2,...
    'plineloc',2,...
    'mlabelparallel','north',...
    'plabelmeridian','west',...
    'parallellab','on',...
    'meridianlab','on',...
    'mlabelround',0,...
    'plabelround',0,...
    'grid','on',...
    'flinewidth',1,...
    'frame','on',...
    'labelformat','none');
colormap(jet)
set(gca,'visible','off',...
    'clim',pvals{pidx,1},...
    'xlim',xl,...
    'ylim',yl)
lx=lonlim(1)+0.02*diff(lonlim);
ly=latlim(1)+0.05*diff(latlim);

[xr,yr]=projfwd(getm(gca),ly,lx);
rh = scaleruler;
setm(rh,'units','km',...
    'majortick',(0:50:150),...
    'majorticklength',5,...
    'minortick',[0 0],...
    'minorticklength',0,...
    'xloc',xr,...
    'yloc',yr,...
    'color','k',...
    'linewi',2,...
    'fontang','it',...
    'fontweight','b')
    
c2=colorbar('horiz');
cp=get(c2,'position');
cp2=cp;
cp2(1)=cp(1)+0.425;
cp2(2)=cp(2)+0.7;
cp2(3)=cp(3)*0.3;
cp2(4)=cp(4)*0.8;
set(c2,'position',cp2);
set(get(c2,'xlabel'),'string',...
    pvals{pidx,2})
