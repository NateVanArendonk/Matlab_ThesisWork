% Load NRRP
% Variables of Interest are:
%   PSFC - Surface pressure [Pa]
%   U10 - Wind at 10m, U dir [m/s]
%   V10
%   MAXWSPD - max wndspd at 10m [m/s]
%
% Note that Pressure, U, and V are ALL on separate grids. (Arakawa C-grid)
%   Press is on the "Mass" grid
%   U on U grid
%   V on V grid
%   position infor in 'geo_em.d02.nc' (nested I believe)
clearvars
return
% Define sites of interest to extract
%fname = 'D:\GoogleDrive\CoastalChangeModeling\Workplan\WeatherStationMeta.csv';
%Station = readStationMeta( fname );

%sites = [1:12 14:23];

% Load Model domain grids
fol_loc = 'C:\Users\ahooshmand\Desktop\Data_Forcings\Wind_Products\ExtractNNRP\';
domaininfo = ncinfo('geo_em.d02.nc');
lat_m = ncread('geo_em.d02.nc','XLAT_M');
lon_m = ncread('geo_em.d02.nc','XLONG_M');
lat_u = ncread('geo_em.d02.nc','XLAT_U');
lon_u = ncread('geo_em.d02.nc','XLONG_U');
lat_v = ncread('geo_em.d02.nc','XLAT_V');
lon_v = ncread('geo_em.d02.nc','XLONG_V');
lat_c = ncread('geo_em.d02.nc','CLAT');     %Computational grid (same as mass grid)
lon_c = ncread('geo_em.d02.nc','CLONG');
elev_m = ncread('geo_em.d02.nc','HGT_M');    %Height of topograpghy in meters MSL
land_mask = ncread('geo_em.d02.nc','LANDMASK');
[N, M] = size(lat_m);
return
% Find nearest output grid cells
min_dist = 15000; %[m]
for nn = 1:length(sites)
    pos_lat = Station.lat(sites(nn));
    pos_lon = Station.lon(sites(nn));
    [ i_loc(nn), j_loc(nn) ] = findNearestGridPos( lat_m, lon_m, pos_lat, pos_lon, min_dist );
end

% Load coarse bathy
Bathy = ncdataset('etopo1_bedrock_PugetSound.nc');
bathy_lat = Bathy.data('lat');
bathy_lon = Bathy.data('lon');
bathy_Band1 = Bathy.data('Band1');

% Set times to load
time_vec = datenum(1950,1,1):(1/4):datenum(2010,12,31,18,0,0);

% Preallocate
slp = NaN(length(time_vec),length(sites));
wndspd = slp;
wnddir = slp;
sst = slp;
surfp = slp;
maxwspd = slp;

tic
count = 0;
for tt=1:length(time_vec)
    
    % Load model prediction time step
    fol_loc = sprintf('../NNRP/%d/',year(time_vec(tt)));
    file_name = sprintf('wrfoutp_d02_%04d%02d%02d%02d.nc',...
        year(time_vec(tt)),month(time_vec(tt)),day(time_vec(tt)),hour(time_vec(tt)));
    model_path = [fol_loc file_name];
    modelinfo = ncinfo(model_path);
    
    % Load in desired variables
    time = ncread(model_path,'Time');
    surf_press = ncread(model_path,'PSFC');
    u10 = ncread(model_path,'U10');
    v10 = ncread(model_path,'V10');
    maxwspd_mat = ncread(model_path,'MAXWSPD');
    %skinTemp = ncread(model_path,'SKINTEMP'); %[K]
    sst_mat = ncread(model_path,'SST')+273.15; %[C]
    Temp2m = ncread(model_path,'T2'); %[K]
    
    % Convert surface pressure to SLP
    surf_press_adj = convertToSLP(surf_press,elev_m,Temp2m);
    
    % Extract just neareast grid cells
    slp(tt,:) = double(surf_press_adj(sub2ind([N, M],i_loc,j_loc)));
    surfp(tt,:) = double(surf_press(sub2ind([N, M],i_loc,j_loc)));
    u10 = double(u10(sub2ind([N, M],i_loc,j_loc)));
    v10 = double(v10(sub2ind([N, M],i_loc,j_loc)));
    sst(tt,:) = double(sst_mat(sub2ind([N, M],i_loc,j_loc)));
    maxwspd(tt,:) = double(maxwspd_mat(sub2ind([N, M],i_loc,j_loc)));
    
    % Calc wind speed
    wndspd(tt,:) = sqrt(u10.^2 + v10.^2);
    wnddir_temp = (180/pi)*atan2(v10,u10);
    
    % Rotate winds to compass directions
    wnddir_temp = 90 - wnddir_temp;
    wnddir_temp(wnddir_temp<0)=wnddir_temp(wnddir_temp<0)+360;
    
    % switch to conventional coming from dir
    wnddir_temp = wnddir_temp+180;
    wnddir(tt,:) = wrapTo360(wnddir_temp);
    
    count = count+1;
    if count > 1000
        count = 0;
        toc
        tic
    end
    
end
toc


% Save
for nn = 1:length(sites)
    O.slp = slp(:,nn);
    O.surfp = surfp(:,nn);
    O.sst = sst(:,nn);
    O.maxwspd = maxwspd(:,nn);
    O.wndspd_10m = wndspd(:,nn);
    O.wnddir = wnddir(:,nn);
    O.time = time_vec(:,nn);
    outname = sprintf('Output/NNRP_%s.mat',Station.shortName{sites(nn)});
    save(outname,'-struct','O');
end

return



% Average U, V gridded winds to Mass grid (ALREADY AVERAGED)
% u_avg = 0.5*diag(ones(N+1,1))+0.5*diag(ones(N,1),1);
% u_avg(:,end) = [];
% u10 = u_avg*u10;

% % Read in just one grid point (may not be faster)
% tic
% surf_press1 = ncread(model_path,'PSFC',[60 60 1],[1 1 1]);
% toc

%% Initial plotting (over correcting pressure)
clf
subplot(121)
pcolor(lon_m,lat_m,surf_press)
shading flat
hold on
caxis([97000 105000])
xlim([-124 -122])
ylim([47 49])
contour(bathy_lon,bathy_lat,bathy_Band1,[0 0],'k')
colorbar

subplot(122)
pcolor(lon_m,lat_m,surf_press)
shading flat
hold on
caxis([97000 105000])
xlim([-124 -122])
ylim([47 49])
contour(bathy_lon,bathy_lat,bathy_Band1,[0 0],'k')
colorbar
colormap(jet)
% Variables of interest
% U10, V10, MAXWSPD, Time





%% FIND Best grid output locatoins (COAST

sites = 1:23;
clf
inds = land_mask == 0;
plot(lon_m,lat_m,'o','Color',[.6 .6 .6])
hold on
plot(lon_u,lat_u,'>','Color',[.6 .6 .6])
plot(lon_v,lat_v,'^','Color',[.6 .6 .6])
%plot(lon_c,lat_c,'*','Color',[.6 .6 .6])
plot(lon_m(inds),lat_m(inds),'o','Color',lines(1),'MarkerFaceColor',lines(1))
%plot(lon_u(inds),lat_u(inds),'>','Color',lines(1),'MarkerFaceColor',lines(1))
%plot(lon_v(inds),lat_v(inds),'^','Color',lines(1),'MarkerFaceColor',lines(1))


shading flat
caxis([95000 103000])
xlim([-125 -122])
ylim([46 49])
contour(bathy_lon,bathy_lat,bathy_Band1,[0 0],'k')
plot(Station.lon(sites),Station.lat(sites),'m*','MarkerSize',12)
colormap(bone)
for nn=sites
    text(Station.lon(nn),Station.lat(nn),Station.shortName{nn})
end
plot(lon_m(sub2ind(size(lat_m),i_loc,j_loc)),lat_m(sub2ind(size(lat_m),i_loc,j_loc)),'or','MarkerSize',12)

fout_name = 'NNRP_grids_WA_Coast';
hFig = gcf;
hFig.PaperUnits = 'inches';
hFig.PaperSize = [8.5 11];
hFig.PaperPosition = [0 0 8.5 11];
print(hFig,'-dpdf',fout_name)


%% FIND Best grid output locatoins ZOOOM

sites = 1:23;
clf
inds = land_mask == 0;
plot(lon_m,lat_m,'o','Color',[.6 .6 .6])
hold on
plot(lon_u,lat_u,'>','Color',[.6 .6 .6])
plot(lon_v,lat_v,'^','Color',[.6 .6 .6])
%plot(lon_c,lat_c,'*','Color',[.6 .6 .6])
plot(lon_m(inds),lat_m(inds),'o','Color',lines(1),'MarkerFaceColor',lines(1))
%plot(lon_u(inds),lat_u(inds),'>','Color',lines(1),'MarkerFaceColor',lines(1))
%plot(lon_v(inds),lat_v(inds),'^','Color',lines(1),'MarkerFaceColor',lines(1))


shading flat
caxis([95000 103000])
xlim([-124 -122])
ylim([47 49])
contour(bathy_lon,bathy_lat,bathy_Band1,[0 0],'k')
plot(Station.lon(sites),Station.lat(sites),'m*','MarkerSize',12)
colormap(bone)
for nn=sites
    text(Station.lon(nn),Station.lat(nn),Station.shortName{nn})
end
plot(lon_m(sub2ind(size(lat_m),i_loc,j_loc)),lat_m(sub2ind(size(lat_m),i_loc,j_loc)),'or','MarkerSize',12)

fout_name = 'NNRP_grids_ZOOM';
hFig = gcf;
hFig.PaperUnits = 'inches';
hFig.PaperSize = [8.5 11];
hFig.PaperPosition = [0 0 8.5 11];
print(hFig,'-dpdf',fout_name)

