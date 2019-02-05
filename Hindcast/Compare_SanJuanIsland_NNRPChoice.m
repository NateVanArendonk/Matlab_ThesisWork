% Load Coastline
load('WA_coast.mat');
% Find all shapes above threshold
wa_lat = [];
wa_lon = [];
thresh = 500 ;
for j = 1:length(wa_coast)
    temp_x = wa_coast(j).X;
    temp_y = wa_coast(j).Y;
    if length(temp_x) >= thresh %&& j ~= 3348 % 3348 is oregon
        for m = 1:length(temp_x);
            wa_lat(end+1) = temp_y(m);
            wa_lon(end+1) = temp_x(m);
        end
    end
end

% NNRP Points
np = dir('E:\Abbas\PS_COSMOS\Thesis_Modeling\Quantile_Correction\NNRP_WaterPointData\*.mat');
for f = 1:length(np)
    f_open = strcat('E:\Abbas\PS_COSMOS\Thesis_Modeling\Quantile_Correction\NNRP_WaterPointData\',np(f).name);
    n = load(f_open);
%     nn(f).lat = n.Position(2);
%     nn(f).lon = n.Position(1);
    nn(f).lat = n.lat;
    nn(f).lon = n.lon(1);
    nn(f).name = np(f).name;
end
return

%%

%% Find Water Points
water_points = [water_lon,water_lat];
waterInds = [];
for ii = 1:length(nn)
    cur_points = [nn(ii).lon,nn(ii).lat];
    in = ismember(cur_points,water_points);
    if sum(in) == 2
        waterInds(end+1) = ii;
    end
end
    
    
nn_water = nn(waterInds);

%% Copy water points to water point file
file_names = cell(length(np),1);
for ii = 1:length(file_names)
    file_names{ii} = np(ii).name;
end
for ii = 1:length(nn_water)
    for jj = 1:length(file_names)
        if strcmp(nn_water(ii).name,file_names{jj})
            filename = file_names{jj};
            copyfile(['NNRP_PointData\' filename], 'NNRP_WaterPointData')
        end
    end
end
            
    