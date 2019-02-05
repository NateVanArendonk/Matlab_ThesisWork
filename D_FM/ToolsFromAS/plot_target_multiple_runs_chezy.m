base_dir='y:\Projects\pscosmos\runs\roughness_sensitivity\';
mruns=cellfun(@(x)([sprintf('psfm_r%0.0f',x),'\DFM_OUTPUT_ps_2d\']),...
    num2cell(([4 2 1 3 5])'),'un',0);
fname='ps_2d_his.nc';


info=cellfun(@(x)(ncinfo([base_dir,x,fname])),mruns);
vars=arrayfun(@(y)(arrayfun(@(x)(x.Name),y.Variables,'un',0)'),...
    info,'un',0);
stns=cellfun(@(x)(cellstr(ncread([base_dir,x,fname],'station_name')')),...
    mruns,'un',0);

%grab model data for specified stations
dstns={'NOAA_9443090_Neah_Bay';...
    'NOAA_9444090_Port_Angeles';...
    'NOAA_9444900_Port_Townsend';...
    'NOAA_9447130_Seattle';...
    'NOAA_9446484_Tacoma';...
    'NOAA_9449424_Cherry_Point'};
[~,~,istns]=cellfun(@(x)(intersect(dstns,x,'stable')),stns,'un',0);
gstns=cellfun(@(x,y)(x(y)),stns,istns,'un',0);


tstart=datenum(2017,2,1);
tend=datenum(2017,4,1);

time=cellfun(@(x)(ncread([base_dir,x,fname],'time')),mruns,'un',0);

it=datenum(2017,1,1); %from netcdf attributes
mtime=cellfun(@(x)(it+(x./86400)),time,'un',0);

st=cellfun(@(x)(find(x>=tstart,1,'first')),mtime,'un',0);
ei=cellfun(@(x)(find(x<=tend,1,'last')),mtime,'un',0);
mtime=cellfun(@(x,y,z)(x(y:z)),mtime,st,ei,'un',0);

wl=cellfun(@(z,s,n,i)(cellfun(@(x)(ncread([base_dir,z,fname],'waterlevel',...
    [x s],[1 numel(n)])'),num2cell(i),'un',0)),mruns,st,mtime,istns,'un',0);
lon=cellfun(@(z,s,n,i)(cellfun(@(x)(ncread([base_dir,z,fname],...
    'station_x_coordinate',[x 1],[1 1])'),num2cell(i),'un',0)),...
    mruns,st,mtime,istns,'un',0);
lat=cellfun(@(z,s,n,i)(cellfun(@(x)(ncread([base_dir,z,fname],...
    'station_y_coordinate',[x 1],[1 1])'),num2cell(i),'un',0)),...
    mruns,st,mtime,istns,'un',0);


dname='y:\Projects\pscosmos\data\wl\coops_2017\navd88\';
fnames=dir([dname,'*.mat']);
corder=[1 2 3 5 4 6];
wld=arrayfun(@(x)(x.wld),...
    arrayfun(@(x)(load([dname,x.name])),fnames(corder)));
wl_raw=arrayfun(@(x)(double(x.WL_VALUE)),wld,'un',0);
[wld(:).WL_VALUE]=deal(wl_raw{:});

%interp data onto model timeframe
wli=arrayfun(@(x)(interp1(x.time,x.WL_VALUE,mtime{1})),wld,'un',0);

%calculate target stats
dmwl=cellfun(@mean,wli);

mwl=cellfun(@(x)(cell2mat(x')),wl,'un',0)';
mmwl=cellfun(@(x)(mean(x,1)'),mwl,'un',0);

bias=cellfun(@(x)(x-dmwl),mmwl,'un',0)';

dwlu=cellfun(@(x,y)(x-y),wli,num2cell(dmwl),'un',0); %remove mean from obs
crmsd=cell(length(bias),1);
rmsd=cell(length(bias),1);
for i=1:length(bias)
    mwlu=num2cell(bsxfun(@minus,mwl{i},mmwl{i}'),1)'; %remove mean from mod
    crmsd{i}=cellfun(@(x,y)(sqrt((sum((x-y).^2))/numel(x))),mwlu,dwlu).*... %craziness
        cellfun(@(x,y)(sign(std(x)-std(y))),num2cell(mwl{i},1)',wli); 
    rmsd{i}=cellfun(@(x,y)(sqrt((sum((x-y).^2))/numel(x))),...
        num2cell(mwl{i},1)',wli)';
end


%plot the target diagram
crmsdm=cell2mat(crmsd');
biasm=cell2mat(bias');
maxr=max([crmsdm(:);biasm(:)]);
maxr2=maxr+0.05*maxr;

xi=linspace(0,2*pi,100);
r=0:0.05:0.25;
[x,y]=cellfun(@(x)(pol2cart(xi,x)),num2cell(r),'un',0);



f=figure;
set(f,'renderer','zbuffer','units','inches',...
    'position',[1 1 8 7],...
    'paperpositionmode','au',...
    'color','w','inverthard','off')

hold on 
cellfun(@(x,y)(plot(x,y,'k--')),x,y)

plot(crmsdm',biasm','color',[0.6 0.6 0.6])

% th=text(crmsdmat(1,:),bmat(1,:),stns,...
%     'fontang','it','fontweight','b','horizontalalign','right',...
%     'verticalalign','top')
% set(th([1,3]),'horizontalalign','left')



mkrs={'o';'^';'s';'>';'p';'<'};
for i=1:length(mkrs)
    dh(i)=plot(999,999,'marker',mkrs{i},...
        'color','k','markerfacecolor','k',...
        'linestyle','none');
end

cols=num2cell(jet(length(mruns)),2);
for i=1:length(dstns)
    for j=1:length(mruns)
   lh(i,j)=plot(crmsdm(i,j),biasm(i,j),...
        'color',cols{j},'marker',mkrs{i},...
        'markerfacecolor',cols{j},'linestyle','none');
    end

end


lh2=lh(1,:);
lht=[dh(:);lh2(:)];

legstr={'Chezy = 55';...
    'Chezy = 60'
    'Chezy = 62.651';...
    'Chezy = 65';...
    'Chezy = 70'};
leg=legend(lht,cat(1,gstns{1},legstr));
set(leg,'location','northeastoutside',...
    'fontang','it','fontweight','b',...
    'box','off','interpreter','none')

set(gca,'da',[1 1 1],...
    'xaxislocation','origin',...
    'yaxislocation','origin',...
    'xlim',[-0.3 0.3],...
    'ylim',[-0.3 0.3],...
    'xtick',r,'ytick',r)

xl=xlim;
text(xl(2),0,'\bf\itRMSD''(m)',...
    'horizontalalign','right',...
    'verticalalign','bot')
text(0,xl(2),'\bf\itBias (m)',...
    'horizontalalign','right',...
    'rotation',90,...
    'verticalalign','top')








