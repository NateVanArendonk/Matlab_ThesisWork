dname='y:\Projects\pscosmos\data\salinity\';
fnames=dir([dname,'*2015.csv']);

lfile=[dname,'stn_locs.txt'];
fid=fopen(lfile,'r');
sdata=textscan(fid,'%s%f%f','delimiter','\t');
fclose(fid);

dstr=repmat(struct('station',[],...
    'lon',[],...
    'lat',[],...
    'mtime',[],...
    'depth',[],...
    'temp',[],...
    'sal',[]),length(fnames),1);

for i=1:length(fnames)
    fmt='%s%s%*f%f%f%f%*[^\n]';
    fid=fopen([dname,fnames(i).name],'r');
    data=textscan(fid,fmt,'headerlines',1,...
        'delimiter',',');
    fclose(fid);
    
    dn=datenum(data{2});
    udates=unique(dn);
    fun=@(x)([min(x) max(x)]);
    dz=0.5;
    
    zlim=fun(data{3});
    z=(zlim(1):dz:zlim(2))';
    zm=repmat(z,1,numel(udates));
    tm=repmat(udates',numel(z),1);
    
    s=scatteredInterpolant(dn,data{3},data{4});
    salm=s(tm,zm);
    t=scatteredInterpolant(dn,data{3},data{5});
    tempm=t(tm,zm);
    
%     rho = gsw_rho(salm,tempm,zm);
    
    stn=data{1}{1}(2:end-1);
    sidx=strmatch(stn,sdata{1});
    dstr(i).station=sdata{1}(sidx);
    dstr(i).lon=sdata{2}(sidx);
    dstr(i).lat=sdata{3}(sidx);
    dstr(i).depth=-z;
    dstr(i).mtime=udates;
    dstr(i).sal=salm;
    dstr(i).temp=tempm;
%     dstr(i).rho=rho;
end

msal=cellfun(@mean,arrayfun(@(x)(mean(x.sal)),dstr,'un',0));
mtemp=cellfun(@mean,arrayfun(@(x)(mean(x.temp)),dstr,'un',0));
lon=[dstr(:).lon]';
lat=[dstr(:).lat]';
dlmwrite('c:\data\projects\pscosmos\salinity\ps2015_sal_avg.xyz',...
    [lon,lat,msal],'delimiter','\t','precision','%0.7e')



