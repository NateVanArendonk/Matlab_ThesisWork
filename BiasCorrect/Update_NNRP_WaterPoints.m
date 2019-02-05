% Load in Old NNRP Water Points and get list of all Lat Lon points 
clearvars 
% Get list of NNRP Points
d = dir('E:\Abbas\PS_COSMOS\Thesis_Modeling\Quantile_Correction\NNRP_WaterPointData\*.mat');
fol = 'E:\Abbas\PS_COSMOS\Thesis_Modeling\Quantile_Correction\NNRP_WaterPointData\';

% Loop through and save Lat Lon Points in to a structure 
for ii = 1:length(d)
    temp = load([fol d(ii).name]);
    N(ii).name = d(ii).name;
    N(ii).lon = temp.lon;
    N(ii).lat = temp.lat;
end
% save('NNRP_WaterPointsStruct','N');
%% Load in West Coast NNRP matfile 
M = load('west_coast_NNRP.mat');
%%
for f = 1:length(N)
    [r,c] = find(N(f).lat == M.lat_m & N(f).lon == M.lon_m);
    lat = N(f).lat;
    lon = N(f).lon;
    % Grab data at exact location
    elev = M.elev_m(r,c);
    slp = squeeze(M.slp(r,c,:));
    wndspd = squeeze(M.wndspd(r,c,:));
    wnddir = squeeze(M.wnddir(r,c,:));
    surfp = squeeze(M.surfp(r,c,:));
    wind_u = squeeze(M.wind_u(r,c,:));
    wind_v = squeeze(M.wind_v(r,c,:));
    time = M.time;
    
    % Save the data
    outname = N(f).name;
    save(outname,'elev','slp','wndspd','wnddir','surfp','time','lat','lon','wind_u','wind_v')
    movefile(outname,'NNRP_Updated\')
    
    fprintf('Completed %d out of %d\n',f,length(N))
end