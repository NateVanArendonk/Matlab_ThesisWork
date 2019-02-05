% Extract NNRP Points from NNRP Point List

% Load in Data for Quantile Correction 
% station_name = 'tacoma_narrows_arpt_hourly';
% data_fol = 'C:\Users\ahooshmand\Desktop\Data_Forcings\station_data\gap_hourly'; % Folder of Obs Data
% O = load(strcat(data_fol,'\',station_name));
% N = load('C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\Storm_Tracker\west_coast_NNRP.mat');
% n_pt = load('NNRP_Points\South_Tacoma_NNRP_Point');


clearvars

N = load('WA_NNRP.mat'); % Load in West Coast NNRP 
%%
np = dir('nn2save\*.mat'); % List of NNRP Points that I want to save
for f = 1:length(np)
    fol_name = 'NNRP_PointData\';
    if ~exist([fol_name np(f).name])
        fprintf('Saving %s\n',np(f).name)
        f_open = strcat('nn2save\',np(f).name);
        [~,name,~] = fileparts(f_open);
        n = load(f_open);
        % Find position on NNRP Grid
        if isfield(n,'plat')
            [r,c] = find(N.lat_m == n.plat & N.lon_m == n.plon);
            lat = n.plat;
            lon = n.plon;
        else
            [r,c] = find(N.lat_m == n.Position(2) & N.lon_m == n.Position(1));
            lat = n.Position(2);
            lon = n.Position(1);
        end
        % Grab data at exact location
        elev = N.elev_m(r,c);
        slp = squeeze(N.slp(r,c,:));
        wndspd = squeeze(N.wndspd(r,c,:));
        wnddir = squeeze(N.wnddir(r,c,:));
        surfp = squeeze(N.surfp(r,c,:));
        time = N.time;

        
        % Save the data
        outname = sprintf('%s',name);
        save(outname,'elev','slp','wndspd','wnddir','surfp','time','lat','lon')
    end
end

