dname='y:\Projects\pscosmos\runs\psfm_3\DFM_OUTPUT_salishc\';
fname='salishc_his.nc';

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


tstart=datenum(2017,8,3);
tend=datenum(2017,10,1);
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

mstruct=gcm(axesm('tranmerc'));
mstruct.falsenorthing=0;
mstruct.falseeasting=500000;
mstruct.geoid=[6378137 0.081819191042815];
mstruct.mapparallels=0;
mstruct.nparallels=1;
mstruct.origin=[0 -123 0];
mstruct.scalefactor=0.9996;
[x,y]=projfwd(mstruct,cell2mat(lat),cell2mat(lon));
dist=[0;cumsum(sqrt(diff(x).^2+diff(y).^2))];

xc=num2cell(x);
yc=num2cell(y);
dc=num2cell(dist);

wl=cellfun(@(x)(ncread([dname,fname],'waterlevel',...
    [x st],[1 numel(mtime)])'),num2cell(istns),'un',0);


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

%get the data

%grab model data for specified stations
gstns={'NOAA_9443090_Neah_Bay';...
    'NOAA_9444090_Port_Angeles';...
    'NOAA_9444900_Port_Townsend';...
    'NOAA_9447130_Seattle';...
    'NOAA_9446484_Tacoma'};



% %grab data from noaa coops website
% cnames=cellfun(@(x)(x(6:12)),gstns,'un',0);
% opt.start_time=datestr(mtime(1));
% opt.end_time=datestr(mtime(end));
% opt.datum='MSL';
% opt.out_dir='z:\projects\pscosmos\data\';
% wld=cellfun(@(x)(get_coops_erdapp(x,opt)),cnames);

cnames=cellfun(@(x)(x(6:12)),gstns,'un',0);
dval='IOOS_SixMin_Verified_Water_Leve';
dname='z:\projects\pscosmos\data\wl\';
for i=1:length(cnames)
    cfiles=dir([dname,'*',cnames{i},'*.mat']);
    fnames={cfiles(:).name}';
    data=arrayfun(@(x)(x.(dval)),...
        cellfun(@(x)(load([dname,x])),fnames));
    %concatenate files if necessary
    if length(fnames)>1
        fields=fieldnames(data);
        dc=struct2cell(data)';
        dc2=cellfun(@(x)(cell2mat(dc(:,x))),...
            num2cell(1:length(fields)),'un',0);
        datar=cell2struct(dc2',fields);
        
        %make sure times are unique
        [~,tidx]=unique(datar.time);
        wld(i)=structfun(@(x)(x(tidx,:)),datar,'un',0);
    end
    
    %opendap to datenum
    wld(i).time=(wld(i).time./86400) + 719529;
end
wld=wld(:);


%convert coords to cartesian
[xo,yo]=arrayfun(@(x)(projfwd(mstruct,x.latitude(1),...
    x.longitude(1))),wld);
%calculate nearest point 
dstr=ipdm([xo yo],[x y],'subset','nearest',...
    'result','structure');
disto=dist(dstr.columnindex);


%correct wl data to NAVD88 (same as model bathymetry)
%some station pages do not include published offsets
%some found here: https://www.ngs.noaa.gov/Tidal_Elevation/
%No offset found for Friday Harbor, so removed from analysis for now
offset=[1.05;...  %noaa datums page
    1.1650;...    %noaa datums page
    1.262;...     %noaa tidal datum page (pid: AI2202)
    1.3090;...    %noaa datums page
    1.3560];      %noaa tidal datum page (pid: AI2204)

wl_msl=arrayfun(@(x)(double(x.WL_VALUE)),wld,'un',0);
wl_navd=cellfun(@(x,y)(x+y),wl_msl,num2cell(offset),'un',0);
[wld(:).DATUM]=deal('NAVD88');
[wld(:).WL_VALUE]=deal(wl_navd{:});

%interp data onto model timeframe
wli=arrayfun(@(x)(interp1(x.time,x.WL_VALUE,mtime)),wld,'un',0);




%tidal analysis
const={'M2';'K1';'O1';'N2';'S2'};

mt=arrayfun(@(x)(t_tide(x.wl,'start time',x.time(1),...
    'interval',diff(mtime(1:2))*24,...
    'latitude',x.lat)),mstr);

dlat=arrayfun(@(x)(x.latitude(1)),wld,'un',0);
dlon=arrayfun(@(x)(x.longitude(1)),wld,'un',0);
dt=cellfun(@(x,y)(t_tide(x,'start time',mtime(1),...
    'interval',diff(mtime(1:2))*24,...
    'latitude',y)),wli,dlat);




damp=cell(1,length(const));
dpha=cell(1,length(const));
mamp=cell(1,length(const));
mpha=cell(1,length(const));

for j = 1:length(const)
        dtidx=arrayfun(@(x)(find(strcmpi(const{j},...
            cellstr(x.name)))),dt,'un',0);
        [dt(:).midx]=deal(dtidx{:});
        damp{j}=arrayfun(@(x)(x.tidecon(x.midx,1)),dt);
        dpha{j}=arrayfun(@(x)(x.tidecon(x.midx,3)),dt);
        
    mtidx=arrayfun(@(x)(find(strcmpi(const{j},...
        cellstr(x.name)))),mt,'un',0);
    [mt(:).midx]=deal(mtidx{:});
    mamp{j}=arrayfun(@(x)(x.tidecon(x.midx,1)),mt);
    mpha{j}=arrayfun(@(x)(x.tidecon(x.midx,3)),mt);
    
    %fix phases so they don't rotate through 0
    bidx=find(abs(diff(mpha{j}))>300)+1;
    mpha{j}(bidx:end)=mpha{j}(bidx:end)+360;
    
    bidx=find(abs(diff(dpha{j}))>300)+1;
    dpha{j}(bidx:end)=dpha{j}(bidx:end)+360;
    
end





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

cidx=1;
fun=@(x)([min(x) max(x)]);
clims=fun(mamp{cidx});

f=figure;
set(f,'renderer','zbuffer','units','inches',...
    'position',[0 0 9 7],...
    'paperpositionmode','au')

subplot(1,3,1)
axesm mercator
hold on 
geoshow(A,R)
scatterm(cell2mat(lat),...
    cell2mat(lon),10,mamp{cidx},'filled')

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
    ['\bf\it',const{cidx},' Amp. (m)'])


ax2(1)=subplot(2,3,2:3);
plot([mstr(:).dist]./1000,mamp{cidx},'ko-')
hold on 
 plot(disto./1000,damp{cidx},'ro','markerfacecolor','r')
 leg=legend('Model','Measurements');
set(leg,'location','northwest')
set(gca,'xticklabel',[])
      
    ylabel(['\bf\it',const{cidx},' Amp. (m)'])

ax2(2)= subplot(2,3,5:6);
plot([mstr(:).dist]./1000,mpha{cidx},'ko-')
hold on 
 plot(disto./1000,dpha{cidx},'ro','markerfacecolor','r')
ylabel(['\bf\it',const{cidx},' Phase (deg)'])
xlabel('\bf\itDistance (km)')

set(ax2,'yaxislocation','right')
 
 
 ap=get(ax2,'position');
 offset=num2cell([-0.075;0.075]);
 ap2=cellfun(@(x,y)([x(1) x(2)+y x(3) x(4)*0.8]),ap,offset,'un',0);
 set(ax2,{'position'},ap2)
 