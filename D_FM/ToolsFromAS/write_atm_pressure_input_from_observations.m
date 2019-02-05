

base_url='https://dods.ndbc.noaa.gov/thredds/dodsC/data/stdmet/';
url={'46087/46087h2017.nc';...
    '46088/46088h2017.nc';...
    'wpow1/wpow1h2017.nc';...
    'tcnw1/tcnw1h2017.nc'};

tstart=datenum(2017,1,1);
tend=datenum(2017,4,2);
time=cellfun(@(x)(double(ncread([base_url,x],'time'))),url,'un',0);
dmtime=cellfun(@(x)(datenum(1970,1,1)+(x./86400)),time,'un',0);
st=cellfun(@(x)(find(x>=tstart,1,'first')),dmtime,'un',0);
et=cellfun(@(x)(find(x<=tend,1,'last')),dmtime,'un',0);

dctime=cellfun(@(x,y,z)(x(y:z)),dmtime,st,et,'un',0);

dvals={'air_pressure'};
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

xi=tstart:1/24:tend;
xc=xi(2:end)-(diff(xi)/2);
[n,~,bin]=arrayfun(@(x)(histcounts(x.mtime,xi)),dstr,'un',0);
p=arrayfun(@(x)(x.air_pressure),dstr,'un',0)
hm=cellfun(@(x,y,z)(accumarray(x(x~=0),z(x~=0),...
    [numel(y) 1],@mean,NaN)),bin,n,p,'un',0)

hmc=nanmean(cell2mat(hm),2)';

figure
hold on 
arrayfun(@(x)(plot(x.mtime,x.air_pressure)),dstr)
plot(xc,hmc,'k-','linewi',2)
legend(sites)
datetick('x','keeplimits')

it=datenum(2017,1,1);
mins=(xc-it).*1440;

fout='y:\Projects\pscosmos\runs\meteo_ntr_sensitivity\psfm_n5\';
dlmwrite([fout,'ndbc_atm_press.tim'],...
    [mins' hmc'.*100],...
    'delimiter',' ','precision','%0.6f')