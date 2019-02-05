clearvars
dname='E:\DelftFM\AndrewRuns\RoughnessTesting\Manning\r_316\';
fname    = 'ps_2d_map.nc';
grd         = dflowfm.readNet([dname,fname]);
tim = nc_cf_time([dname,fname]);


data  = dflowfm.readMap([dname,fname],numel(tim));
face.mask = true(1,grd.face.FlowElemSize);
[xg,yg] = dflowfm.peri2cell(grd.face.FlowElemCont_x(:,face.mask),...
    grd.face.FlowElemCont_y(:,face.mask));
tri.mask = face.mask(grd.map3);




bname='E:\Abbas\Modeling Resources\PS_DEM\GEBCO\';
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
lonlim=[-130 -122];
latlim=[46.8 51.8];


[lonm,latm]=meshgrid(lon,lat);



% addpath('C:\Program Files\MATLAB\R2017a\toolbox\map\mapproj\')
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
% [IM.xi,IM.yi] = meshgrid(IM.xm,IM.ym);

% test = surf(IM.xi,IM.yi,zeros(size(IM.xi)),IM.im)
% set(test,'edgecolor','none')

hold on

z=0;
p = trisurf(grd.tri(tri.mask,:),xp,yp,zeros(size(xp)),...
    data.face.zwl(grd.map3(tri.mask)));
set(p,'FaceColor','flat','edgecolor','none'); 


clims=[-2 2];

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
    'clim',clims,...
    'xlim',xl,...
    'ylim',yl)
lx=lonlim(1)+0.01*diff(lonlim);
ly=latlim(1)+0.05*diff(latlim);

[xr,yr]=projfwd(getm(gca),ly,lx);
rh = scaleruler;
setm(rh,'units','km',...
    'majortick',(0:50:200),...
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
    'Water Level')
%%
waitforbuttonpress
for i=1:length(tim)
    data  = dflowfm.readMap([dname,fname],i,'zwl',1);
    set(p,'cdata',...
        data.face.zwl(grd.map3(tri.mask)))
    pause(0.1)
end
