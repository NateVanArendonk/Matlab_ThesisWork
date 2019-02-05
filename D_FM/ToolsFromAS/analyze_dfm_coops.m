dname='y:\Projects\pscosmos\runs\roughness_sensitivity\psfm_r1\DFM_OUTPUT_ps_2d\';
fname='ps_2d_his.nc';

info=ncinfo([dname,fname]);
vars=arrayfun(@(x)(x.Name),info.Variables,'un',0)';
stns=cellstr(ncread([dname,fname],'station_name')');

%grab model data for specified stations
dstns={'NOAA_9443090_Neah_Bay';...
    'NOAA_9444090_Port_Angeles';...
    'NOAA_9444900_Port_Townsend';...
    'NOAA_9447130_Seattle';...
    'NOAA_9446484_Tacoma';...
    'NOAA_9449424_Cherry_Point'};
[~,~,istns]=intersect(dstns,stns,'stable');
gstns=stns(istns);

tstart=datenum(2017,2,1);
tend=datenum(2017,4,1);
% spinup=2; %days

time=ncread([dname,fname],'time');
it=datenum(2017,1,1); %from netcdf attributes
mtime=it+(time./86400);

st=find(mtime>=tstart,1,'first');
ei=find(mtime<=tend,1,'last');
mtime=mtime(st:ei);

wl=cellfun(@(x)(ncread([dname,fname],'waterlevel',...
    [x st],[1 numel(mtime)])'),num2cell(istns),'un',0);
lon=cellfun(@(x)(ncread([dname,fname],'station_x_coordinate',...
    [x 1],[1 1])'),num2cell(istns),'un',0);
lat=cellfun(@(x)(ncread([dname,fname],'station_y_coordinate',...
    [x 1],[1 1])'),num2cell(istns),'un',0);

%repackage
mstr=repmat(struct('station',[],'time',mtime,...
    'wl',[]),length(istns),1);
[mstr(:).station]=deal(gstns{:});
[mstr(:).lat]=deal(lat{:});
[mstr(:).lon]=deal(lon{:});
[mstr(:).wl]=deal(wl{:});
% 
% %plot the stations
% serverURL = ['https://basemap.nationalmap.gov/arcgis/services/',...
%     'USGSTopo/MapServer/WMSServer?request=GetCapabilities&service=WMS'];
% info = wmsinfo(serverURL);
% 
% latlim=[46.75 48.9];
% lonlim=[-125.25 -121.7];
% 
% imageLength=1024*2;
% 
% orthoLayer = info.Layer(1);
%  [A, R] = wmsread(orthoLayer, 'Latlim', latlim,...
%             'Lonlim', lonlim, ...
%             'ImageHeight', imageLength,...
%             'ImageWidth', imageLength');

mkrs={'o';'^';'s';'>';'p';'<'};
cols=num2cell(jet(length(mstr)),2);
% 
% f=figure;
% set(f,'renderer','zbuffer','units','inches',...
%     'position',[0 0 8.5 5],...
%     'paperpositionmode','au')
% 
% axesm mercator
% hold on 
% geoshow(A,R)
% 
% lh=cellfun(@(x,y,z,h)(plotm(x,y,'marker',z,'color','k',...
%     'markersize',10,...
%     'markerfacecolor',h)),...
%     lat,lon,mkrs,cols,'un',0);
% set(gca,'visible','off')
% 
% setm(gca,'maplonlimit',lonlim,...
%     'maplatlimit',latlim,...
%     'mlabelloc',1,...
%     'mlineloc',1,...
%     'plabelloc',1,...
%     'plineloc',1,...
%     'mlabelparallel','north',...
%     'plabelmeridian','west',...
%     'parallellab','on',...
%     'meridianlab','on',...
%     'mlabelround',0,...
%     'plabelround',0,...
%     'grid','on',...
%     'flinewidth',2,...
%     'frame','on',...
%     'labelformat','none');
% 
% lx=lonlim(1)+0.01*diff(lonlim);
% ly=latlim(1)+0.05*diff(latlim);
% 
% [xr,yr]=projfwd(getm(gca),ly,lx);
% rh = scaleruler;
% setm(rh,'units','km',...
%     'majortick',(0:25:50),...
%     'majorticklength',5,...
%     'minortick',[0 0],...
%     'minorticklength',0,...
%     'xloc',xr,...
%     'yloc',yr,...
%     'color','k',...
%     'linewi',2,...
%     'fontang','it',...
%     'fontweight','b')
%     
% leg=legend(cat(1,lh{:}),gstns);
% set(leg,'location','northeastoutside',...
%     'fontang','it','fontweight','b',...
%     'box','off','interpreter','none')



%grab data from noaa coops website
% cnames=cellfun(@(x)(x(6:12)),gstns,'un',0);
% opt.start_time=datestr(mtime(1));
% opt.end_time=datestr(mtime(end));
% opt.datum='MSL';
% opt.out_dir='z:\projects\pscosmos\data\wl\';
% wld=cellfun(@(x)(get_coops_erdapp(x,opt)),cnames);

% cnames=cellfun(@(x)(x(6:12)),gstns,'un',0);
% dval='IOOS_SixMin_Verified_Water_Leve';
% dname='y:\Projects\pscosmos\data\wl\';
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
corder=[1 2 3 5 4 6];
wld=arrayfun(@(x)(x.wld),...
    arrayfun(@(x)(load([dname,x.name])),fnames(corder)));
wl_raw=arrayfun(@(x)(double(x.WL_VALUE)),wld,'un',0);
[wld(:).WL_VALUE]=deal(wl_raw{:});
%correct wl data to NAVD88 (same as model bathymetry)
%some station pages do not include published offsets
%some found here: https://www.ngs.noaa.gov/Tidal_Elevation/
%No offset found for Friday Harbor, so removed from analysis for now
% offset=[1.05;...  %noaa datums page
%     1.1650;...    %noaa datums page
%     1.262;...     %noaa tidal datum page (pid: AI2202)
%     1.3090;...    %noaa datums page
%     1.3560;...    %noaa datums page
%     1.3170];      %noaa tidal datum page (pid: AI2204)
% 
% wl_msl=arrayfun(@(x)(double(x.WL_VALUE)),wld,'un',0);
% wl_navd=cellfun(@(x,y)(x+y),wl_msl,num2cell(offset),'un',0);
% [wld(:).DATUM]=deal('NAVD88');
% [wld(:).WL_VALUE]=deal(wl_navd{:});

%interp data onto model timeframe
wli=arrayfun(@(x)(interp1(x.time,x.WL_VALUE,mtime)),wld,'un',0);

%plotting, start with time-series
for i=1:length(wld)
    
    minx=floor(min(mtime));
    maxx=ceil(max(mtime));
    
    xp1=minx+0.025*(maxx-minx);
    xt=14;
    
    clear ax
    f=figure;
    set(f,'renderer','zbuffer','units','inches',...
        'position',[1 1 8.5 5],...
        'color','w','inverthard','off',...
        'paperpositionmode','au')
    ax(1)=subplot(2,3,1:2);
    plot(mtime,mstr(i).wl)
    hold on
    plot(mtime,wli{i},'r')
    
    yl=ylim;
    yp1=yl(1)+0.95*diff(yl);
    text(xp1,yp1,'\bf\it\fontsize{12}A')
    
    leg=legend('Mod.','Obs.');
    set(leg,'orientation','horiz',...
        'box','off','location','southwest',...
        'fontang','it','fontweight','b')
    title(gstns{i},'interpreter','none',...
        'fontang','it','fontweight','b')
    ylabel('\bf\itWater level (m, NAVD88)')
    
    ax2(1)=subplot(2,3,3);
    plot(wli{i},mstr(i).wl,'k.')
    wlmin=min([wli{i};mstr(i).wl]);
    wlmax=max([wli{i};mstr(i).wl]);
    set(gca,'xlim',[wlmin wlmax],...
        'ylim',[wlmin wlmax],...
        'da',[1 1 1],...
        'nextplot','add');
    line([wlmin wlmax],[wlmin wlmax],'color','r',...
        'linewi',2)
    apos=get(ax(1),'position');
    bpos=get(ax2(1),'position');
    
    bpos(2)=apos(2);
    bpos(4)=apos(4);
    set(ax2(1),'position',bpos)
    
    xp2=wlmin+0.025*(wlmax-wlmin);
    yp2=wlmin+0.95*(wlmax-wlmin);
    text(xp2,yp2,'\bf\it\fontsize{12}C')
    
    xlabel('\bf\itObs. WL (m)')
    ylabel('\bf\itMod. WL (m)')
    
    wldiff=mstr(i).wl-wli{i};
    
    ax(2)=subplot(2,3,4:5);
    plot(mtime,wldiff,'k');
    ylabel('\bf\itMod.- Obs. (m)')
    set(ax,'xlim',[minx maxx],...
        'xtick',(minx:xt:maxx))
    datetick('x',2,'keepticks','keeplimits')
    set(ax(1),'xticklabel',[]);
    linkaxes(ax,'x')
    
    yl=ylim;
    yp1=yl(1)+0.95*diff(yl);
    text(xp1,yp1,'\bf\it\fontsize{12}B')
    
    q1=interp1(linspace(0,100,numel(wldiff)),sort(wldiff),1);
    q99=interp1(linspace(0,100,numel(wldiff)),sort(wldiff),99);
    xi=linspace(q1,q99,30);
    [n,bin]=histc(wldiff,xi);
    
    ax2(2)=subplot(2,3,6);
    bh=bar(xi,(n./numel(bin)).*100,'histc');
    set(bh,'facecolor',[0.6 0.6 0.6])
    
    set(ax2(2),'xlim',[q1 q99])
    xlabel(ax2(2),'\bf\itMod. - Obs. (m)')
    ylabel(ax2(2),'\bf\itPercent')
    
    xp3=q1+0.05*(q99-q1);
    yl2=ylim;
    yp3=yl2(1)+0.95*diff(yl2);
    text(xp3,yp3,'\bf\it\fontsize{12}D')
end

%target diagram

dmwl=cellfun(@mean,wli);
mmwl=arrayfun(@(x)(mean(x.wl)),mstr);

bias=mmwl-dmwl;

dwlu=cellfun(@(x)(detrend(x,'constant')),wli,'un',0); %remove mean from obs
mwlu=arrayfun(@(x)(detrend(x.wl,'constant')),mstr,'un',0);

crmsd=cellfun(@(x,y)(sqrt((sum((x-y).^2))/numel(x))).*...
    sign(std(x)-std(y)),mwlu,dwlu);
rmsd=cellfun(@(x,y)(sqrt((sum((x-y).^2))/numel(x))),...
    arrayfun(@(x)(x.wl),mstr,'un',0),wli);

%plot the target diagram
maxr=max(abs([crmsd(:);bias(:)]));
maxr2=ceil((maxr+0.05*maxr).*10)./10;

xi=linspace(0,2*pi,100);
% r=linspace(0,maxr2,5);
r=(0:0.1:0.3);
[x,y]=cellfun(@(x)(pol2cart(xi,x)),num2cell(r),'un',0);


f=figure;
set(f,'renderer','zbuffer','units','inches',...
    'position',[1 1 6 5],...
    'paperpositionmode','au',...
    'color','w','inverthard','off')

hold on 
cellfun(@(x,y)(plot(x,y,'k--')),x,y)


lh=cellfun(@(x,y,z,h)(plot(x,y,'marker',z,'color',h,...
    'markerfacecolor',h)),...
    num2cell(crmsd),num2cell(bias),mkrs,cols,'un',0);

leg=legend(cat(1,lh{:}),gstns);
set(leg,'location','northeastoutside',...
    'fontang','it','fontweight','b',...
    'box','off','interpreter','none')

set(gca,'da',[1 1 1],...
    'xaxislocation','origin',...
    'yaxislocation','origin',...
    'xlim',[-maxr2 maxr2],...
    'ylim',[-maxr2 maxr2])

xl=xlim;
text(xl(2),0,'\bf\itRMSD''(m)',...
    'horizontalalign','right',...
    'verticalalign','bot')
text(0,xl(2),'\bf\itBias (m)',...
    'horizontalalign','right',...
    'rotation',90,...
    'verticalalign','top')



%tidal analysis
const={'M2';'K1';'O1';'N2';'S2'};

mt=arrayfun(@(x)(t_tide(x.wl,'start time',x.time(1),...
    'interval',diff(mtime(1:2))*24,...
    'latitude',x.lat)),mstr);

dlat=arrayfun(@(x)(x.latitude(1)),wld,'un',0);
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


    
end

%visualize constituents



fun=@(x)([min(x)-0.05*(max(x)-min(x)),...
    max(x)+0.05*(max(x)-min(x))]);
for i=1:length(const)
    amax=fun([damp{i};mamp{i}]);
    pmax=fun([dpha{i};mpha{i}]);
    
    f=figure;
    set(f,'renderer','zbuffer','units','inches',...
        'position',[0 0 8 5],...
        'paperpositionmode','au',...
        'color','w','inverthard','off')
    
    ax(1)=subplot(121);
    line(amax,amax,'color','k')
    hold on
    lh=cellfun(@(x,y,z,h)(plot(x,y,'marker',z,'color',h,...
        'markerfacecolor',h)),...
        num2cell(damp{i}),num2cell(mamp{i}),mkrs,cols,'un',0);
    set(ax(1),'da',[1 1 1],...
        'xlim',amax,...
        'ylim',amax)
    
        leg=legend(cat(1,lh{:}),gstns);
        set(leg,'location','northwest',...
            'fontang','it','fontweight','b',...
            'box','off','interpreter','none',...
            'fontsize',6)
        
    xlabel(['\bf\itObs. ',const{i},' Amp. (m)'])
    ylabel(['\bf\itMod. ',const{i},' Amp. (m)'])
    
    ax(2)=subplot(122);
    line(pmax,pmax,'color','k')
    hold on
    cellfun(@(x,y,z,h)(plot(x,y,'marker',z,'color',h,...
        'markerfacecolor',h)),...
        num2cell(dpha{i}),num2cell(mpha{i}),mkrs,cols,'un',0);
        set(ax(2),'da',[1 1 1],...
        'xlim',pmax,...
        'ylim',pmax)
    
    set(ax,'box','on')
    xlabel(['\bf\itObs. ',const{i},' Phase'])
    ylabel(['\bf\itMod. ',const{i},' Phase'])
end




