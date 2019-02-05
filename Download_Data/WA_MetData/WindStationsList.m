stns = dir('C:\Users\ahooshmand\Desktop\Data_Download_New\combined_station_data\*.mat');
for ii = 1:length(stns)
    Z = load(['C:\Users\ahooshmand\Desktop\Data_Download_New\combined_station_data\' stns(ii).name]);
    S(ii).name = Z.name;
    S(ii).lon = Z.lon(1);
    S(ii).lat = Z.lat(1);
end
    