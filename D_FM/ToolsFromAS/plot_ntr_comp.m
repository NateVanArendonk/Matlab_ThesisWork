base_dir='y:\Projects\pscosmos\runs\';
dirs={'psfm_m1\';...
    'psfm_m2\';...
    'psfm_m4\'};
run_id='ps_2d';

fname=cellfun(@(x)([base_dir,x,...
    'DFM_OUTPUT_',run_id,filesep,run_id,'_his.nc']),...
    dirs,'un',0);

info=ncinfo(fname{1}); %assume station list same between runs
vars=arrayfun(@(x)(x.Name),info.Variables,'un',0)';
stns=cellstr(ncread(fname{1},'station_name')');

%grab model data for specified stations
dstns={'NOAA_9443090_Neah_Bay';...
    'NOAA_9444090_Port_Angeles';...
    'NOAA_9444900_Port_Townsend';...
    'NOAA_9447130_Seattle';...
    'NOAA_9446484_Tacoma';...
    'NOAA_9449424_Cherry_Point'};
[~,~,istns]=intersect(dstns,stns,'stable');
gstns=stns(istns);



tstart=datenum(2017,3,3,0,0,0);
tend=datenum(2017,3,31,23,0,0);
% spinup=2; %days

time=cellfun(@(x)(ncread(x,'time')),fname,'un',0);
it=datenum(2017,1,1); %from netcdf attributes
mtime=cellfun(@(x)(it+(x./86400)),time,'un',0);

st=cellfun(@(x)(find(x>=tstart,1,'first')),mtime,'un',0);
ei=cellfun(@(x)(find(x<=tend,1,'last')),mtime,'un',0);
mtime=cellfun(@(x,y,z)(x(y:z)),mtime,st,ei,'un',0);


wl=cellfun(@(y,z,c)(cellfun(@(x)(ncread(y,'waterlevel',...
    [x z],[1 numel(c)])'),num2cell(istns),'un',0)),fname,...
    st,mtime,'un',0);
lon=cellfun(@(x)(ncread([dname,fname],'station_x_coordinate',...
    [x 1],[1 1])'),num2cell(istns),'un',0);
lat=cellfun(@(x)(ncread([dname,fname],'station_y_coordinate',...
    [x 1],[1 1])'),num2cell(istns),'un',0);


cnames=cellfun(@(x)(x(6:12)),gstns,'un',0);

dname='y:\Projects\pscosmos\data\wl\coops_2017\navd88\';
cfiles=dir([dname,'*.mat']);
dnames=arrayfun(@(x)(x.name(6:12)),cfiles,'un',0);
[~,~,ib]=intersect(dnames,cnames)

data=arrayfun(@(x)(x.wld),...
        arrayfun(@(x)(load([dname,x.name])),cfiles(ib)));

    %calculate NTR of data
 addpath('c:\data\tools\CMGmfiles\')
[ntr, jdf]=arrayfun(@(x)(cmglowpass(x.WL_VALUE, 360,...
    x.time,1,'osu')),data,'un',0);
[data(:).ntr]=deal(ntr{:});
[data(:).ntr_dn]=deal(jdf{:});

    
%calculate NTR of model output
[mntr,mntr_dn]=cellfun(@(z,t)(cellfun(@(y)(cmglowpass(y,600,...
    t,1,'osu')),z,'un',0)),wl,mtime,'un',0);


for i=1:length(data)
    figure
    plot(data(i).ntr_dn,data(i).ntr)
    hold on
    cellfun(@(x,y)(plot(x,y)),...
        cellfun(@(x)(x(i)),mntr_dn),cellfun(@(x)(x(i)),mntr));
    
    leg=legend('Observed',...
        'Boundary NTR',...
        'Full Meteo',...
        'No Meteo');
    ylabel('NTR (m, NAVD88)')
    
    set(gca,'xlim',[tstart tend],...
        'xtick',(tstart:7:tend))
    datetick('x',2,'keepticks','keeplimits')
    
    title(gstns(i),'interpreter','none')
end




