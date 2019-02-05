dname='y:\Projects\pscosmos\runs\meteo_ntr_sensitivity\psfm_n2\';
files={'Meteo_2017_int8_slp.nc';...
    'Meteo_2017_int8_u.nc';...
    'Meteo_2017_int8_v.nc'};
info=cellfun(@(x)(ncinfo([dname,x])),files,'un',0);

vars={'slp';...
    'u';...
    'northward_wind'};

lat=ncread([dname,files{1}],'lat');
lon=ncread([dname,files{1}],'lon');

data=cellfun(@(x,y,z)(ncread([dname,x],y,...
    [1 1 1],[size(lat,1) size(lon,2) 1])),...
    files,vars,info,'un',0);
time=ncread([dname,files{1}],'time');
it=datenum(2017,1,1);
mtime=it+(time./24);


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
