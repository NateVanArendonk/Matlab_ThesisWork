% stn='9443090';
% opt.start_time='2016-1-1';
% opt.end_time='2017-1-1';
% opt.datum='NAVD';
% opt.out_dir='z:\projects\pscosmos\data\wl\';
% wld=get_coops_erdapp(stn,opt);

% ddir='y:\Projects\pscosmos\data\wl\coops_2017\raw\neah_bay\';
% dfiles=dir([ddir,'*.mat']);
% fnames={dfiles(:).name}';
% dval='IOOS_SixMin_Verified_Water_Leve';
%     data=arrayfun(@(x)(x.(dval)),...
%         cellfun(@(x)(load([ddir,x])),fnames));
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
%         wld=structfun(@(x)(x(tidx,:)),datar,'un',0);
%     end
%     
%     %opendap to datenum
%     wld.time=(wld.time./86400) + 719529;
%     
%     sidx=6;
%     offset=[1.05;...  %noaa datums page
%     1.1650;...    %noaa datums page
%     1.262;...     %noaa tidal datum page (pid: AI2202)
%     1.3090;...    %noaa datums page
%     1.3560;...    %noaa datums page
%     1.3170];      %noaa tidal datum page (pid: AI2204)
% 
% 
% wld.WL_VALUE=wld.WL_VALUE+offset(sidx);
% wld.DATUM='NAVD88';

load(['y:\Projects\pscosmos\data\wl\coops_2017\navd88\',...
    'NOAA_9443090_Neah_Bay_navd88.mat'])


    
 addpath('c:\data\tools\CMGmfiles\')
[ntr, jdf]=cmglowpass(wld.WL_VALUE, 360, wld.time,1,'osu');

doy=date2doy(jdf);
minx=floor(min(jdf));
maxx=ceil(max(jdf));


figure
plot(jdf,ntr)
hold on 
line([minx maxx],[1.059 1.059],'color','r')
set(gca,'xlim',[minx maxx])
datetick('x','mmm','keeplimits')
ylabel('\bf\itWater level (m, NAVD88)')
xlabel('\bf\it2017')

leg=legend('Non-tidal residual','Station MSL');

wl_corr=ntr;
it=datenum(2017,1,1);
mins=(jdf-it).*1440;

dout='y:\Projects\pscosmos\runs\psfm_6\';
fout='water_level_correction.tim';
dlmwrite([dout,fout],[mins wl_corr-0.069],'delimiter',' ',...
    'precision','%0.7e')


