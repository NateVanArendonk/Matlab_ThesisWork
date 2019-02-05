


base_url='https://dods.ndbc.noaa.gov/thredds/dodsC/data/stdmet/';
url={'46087/46087h2017.nc';...
    '46088/46088h2017.nc';...
    'wpow1/wpow1h2017.nc'};

tstart=datenum(2017,2,1);
tend=datenum(2017,4,1);
time=cellfun(@(x)(double(ncread([base_url,x],'time'))),url,'un',0);
dmtime=cellfun(@(x)(datenum(1970,1,1)+(x./86400)),time,'un',0);
st=cellfun(@(x)(find(x>=tstart,1,'first')),dmtime,'un',0);
et=cellfun(@(x)(find(x<=tend,1,'last')),dmtime,'un',0);

dctime=cellfun(@(x,y,z)(x(y:z)),dmtime,st,et,'un',0);

dvals={'wind_spd';...
    'wind_dir';...
    'air_pressure'};
bdata=cellfun(@(z,s1,s2,n)(cellfun(@(x)(squeeze(ncread(...
    [base_url,z],x,[1 1 s1],...
    [1 1 numel(n)]))),dvals,'un',0)),url,st,et,dctime,'un',0);
bdatac=cat(1,bdata{:});
blat=cellfun(@(x)(ncread([base_url,x],'latitude')),url,'un',0);
blon=cellfun(@(x)(ncread([base_url,x],'longitude')),url,'un',0);
sites=cellfun(@(x)(strtok(x,'/')),url,'un',0);

for i=1:length(url)
    dstr(i).station_name=sites{i};
    dstr(i).longitude=blon{i};
    dstr(i).latitude=blat{i};
    dstr(i).mtime=dctime{i};
    for j=1:length(dvals)
        dstr(i).(dvals{j})=bdata{i}{j};
    end
end

dname='y:\Projects\pscosmos\data\meteo\';
files={'Meteo_2017JFMA_single_slp.nc';...
    'Meteo_2017JFMA_single_u.nc';...
    'Meteo_2017JFMA_single_v.nc'};
info=cellfun(@(x)(ncinfo([dname,x])),files,'un',0);

vars={'slp';...
    'u';...
    'v'};

lat=ncread([dname,files{1}],'lat');
lon=ncread([dname,files{1}],'lon');

dist=arrayfun(@(x)(ipdm([x.longitude x.latitude],[lon(:) lat(:)],...
    'subset','nearest','result','structure')),dstr);
[dx,dy]=arrayfun(@(x)(ind2sub(size(lat),x.columnindex)),dist,'un',0);
dlon=cellfun(@(x,y)(lon(x,y)),dx,dy,'un',0);
dlat=cellfun(@(x,y)(lat(x,y)),dx,dy,'un',0);

time=ncread([dname,files{1}],'time');
it=datenum(2017,1,1);
mtime=it+(time./24);
stm=find(mtime>=tstart,1,'first');
etm=find(mtime<=tend,1,'last');
mctime=mtime(stm:etm);

data=cellfun(@(xi,yi)(cellfun(@(x,y)(squeeze(ncread([dname,x],y,...
    [xi yi 1],[1 1 numel(mctime)]))),...
    files,vars,'un',0)),dx,dy,'un',0);
[mdir,mspd]=cellfun(@(x)(cart2pol(x{2},x{3})),data,'un',0);

figure
subplot(211)
plot(dstr(3).mtime,dstr(3).air_pressure.*100) %data are in hPa
hold on 
plot(mctime,data{3}{1})

subplot(212)
plot(dstr(3).mtime,dstr(3).wind_spd) 
hold on 
plot(mctime,mspd{3})


subplot(313)
plot(dstr(1).mtime,vd{1}) %data are in hPa
hold on 
plot(mctime,data{1}{3})




sf=0.01;
thin=5;

figure
sh=surf(lon,lat,zeros(size(data{1})),data{1});
shading flat
view(2)
hold on 
qh=quiver(lon(1:thin:end,1:thin:end),...
    lat(1:thin:end,1:thin:end),...
    data{2}(1:thin:end,1:thin:end).*sf,...
    data{3}(1:thin:end,1:thin:end).*sf,0);
set(qh,'color','k')

th=title(datestr(mtime(1)));
colorbar

waitforbuttonpress
for i=1:numel(mtime)

    data=cellfun(@(x,y,z)(ncread([dname,x],y,...
    [1 1 i],[size(lat,1) size(lon,2) 1])),...
    files,vars,info,'un',0);

set(sh,'cdata',data{1})
set(qh,'udata',data{2}(1:thin:end,1:thin:end).*sf,...
    'vdata',data{3}(1:thin:end,1:thin:end).*sf);
set(th,'string',datestr(mtime(i)))
pause(0.1)
end
