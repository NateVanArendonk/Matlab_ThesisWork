% ------------------- Load Transects --------------------------------------
kml_fol = 'E:\Abbas\WCRP\Tacoma\Tier3\OrthogTran\Transects\';
k = dir('E:\Abbas\WCRP\Tacoma\Tier3\OrthogTran\Transects\*.kml');
% Load Owen Beach
O = kml2struct([kml_fol k(1).name]);
% Change name of O to have most northly point be 1
num = 0;
count = length(O);
for ii = 1:length(O)
    O(ii).Name = strcat('Line ',' ',num2str(count-num));
    num = num + 1;
end

% Load Ruston and reverse names
R = kml2struct([kml_fol k(2).name]);
lo = length(O);
count = length(R) + lo;
for ii = 1:length(R)
    R(ii).Name = strcat('Line ',' ',num2str(count));
    count = count - 1;
end

% Concatenate Structures
temp = [R,O];
for ii = 1:length(temp)
    T(ii).name = temp(ii).Name;
    T(ii).lon = temp(ii).Lon;
    T(ii).lat = temp(ii).Lat;
    [T(ii).x,T(ii).y] = deg2utm(T(ii).lat,T(ii).lon);
end
clear kml_fol k O num count ii R lo temp

fid = fopen('RW_OB_Transects.kml','w');
fprintf(fid,'<?xml version="1.0" encoding="utf-8"?>\n');
fprintf(fid,'<kml xmlns="http://www.opengis.net/kml/2.2">\n');
fprintf(fid,'   <Document>\n');
fprintf(fid,'      <name>Transect_lines</name>\n');
line_num = length(T);
for ii = 1:length(T)
    fprintf(fid,'      <Placemark>\n');
    fprintf(fid,'         <Snippet maxLines="0"> </Snippet>\n');
    fprintf(fid,'         <description> </description>\n');
    fprintf(fid,'         <name>Line %d</name>\n',line_num);
    fprintf(fid,'         <LineString>\n');
    fprintf(fid,'            <coordinates> %f,%f,0 %f,%f,0</coordinates>\n',T(ii).lon(1),T(ii).lat(1),T(ii).lon(2),T(ii).lat(2));
    fprintf(fid,'         </LineString>\n');
    fprintf(fid,'      </Placemark>\n');
    line_num = line_num - 1;
end
fprintf(fid,   '</Document>\n');
fprintf(fid,'</kml>\n')
fclose(fid);
fclose('all')



%%
% %% OLD UNSURE OF THIS BELOW
% path = 'E:\Abbas\WCRP\Tacoma\Tier3\OrthogTran\Transects\';
%
% dirNames = strcat(path);
% xbr = dir(dirNames); % get list of all runs
% xbr(1:2) = []; % Gets rid of 2 weird results at beginning
% xbr = nestedSortStruct(xbr,{'date'});
% xbr(1) = [];
%
% % ---------- Location of KKL - Load and Convert to UTM
% num_kml = '2';
% switch num_kml % If there are multiple kml files that need to be loaded in
%     case '1'
%         kml_fol = 'E:\Abbas\WCRP\Tacoma\Tier3\OrthogTran\Transects\';
%         kml_nm = 'OwenBeach_Orthog_Transects.kml';
%         temp = kml2struct([kml_fol kml_nm]); % load in KML of Transect
%         for ii = 1:length(temp)
%             T(ii).lat = temp(ii).Lat;
%             T(ii).lon = temp(ii).Lon;
%             [T(ii).x,T(ii).y] = deg2utm(T(ii).lat,T(ii).lon); % Convert to utm
%             T(ii).name = temp(ii).Name;
%         end
%         clear kml_fol kmls kml_nm
%     case '2'
%         % ------------------- Load Transects --------------------------------------
%         kml_fol = 'E:\Abbas\WCRP\Tacoma\Tier3\OrthogTran\Transects\';
%         k = dir('E:\Abbas\WCRP\Tacoma\Tier3\OrthogTran\Transects\*.kml');
%         % Load Owen Beach
%         O = kml2struct([kml_fol k(1).name]);
%         % Change name of O to have most northly point be 1
%         num = 0;
%         count = length(O);
%         for ii = 1:length(O)
%             O(ii).Name = strcat('Line ',' ',num2str(count-num));
%             num = num + 1;
%         end
%
%         % Load Ruston and reverse names
%         R = kml2struct([kml_fol k(2).name]);
%         lo = length(O);
%         count = length(R) + lo;
%         for ii = 1:length(R)
%             R(ii).Name = strcat('Line ',' ',num2str(count));
%             count = count - 1;
%         end
%
%         % Concatenate Structures
%         temp = [R,O];
%         for ii = 1:length(temp)
%             T(ii).name = temp(ii).Name;
%             T(ii).lon = temp(ii).Lon;
%             T(ii).lat = temp(ii).Lat;
%             [T(ii).x,T(ii).y] = deg2utm(T(ii).lat,T(ii).lon);
%         end
%         clear kml_fol k O num count ii R lo temp
% end
%
% % -------------------- Extract Dephts along Transects ---------------------
% % ---------- Name of Folder to house bathy transects
% % transect_folder = 'Owen_Beach_TransectBathy';
%
% % ---------- Create bathy-transects for each transect of KML
% % Set resolution of transect
% ds = 1; %[m]
%
% for ii = 1:length(T)
%     % Create transect
%     [ T(ii).line_x, T(ii).line_y ] = createTransect(T(ii).x(2),T(ii).y(2),T(ii).x(1),T(ii).y(1),ds);
%
%     % Extract subset of depths that covers transect
%     t_inds = dem_x >= min(T(ii).line_x) & dem_x <= max(T(ii).line_x) & ...
%         dem_y >= min(T(ii).line_y) & dem_y <= max(T(ii).line_y);
%     t_x = dem_x(t_inds);
%     t_y = dem_y(t_inds);
%     t_z = dem_z(t_inds);
%
%     % Find nearest depth for each location on transect
%     T(ii).line_z = zeros(size(T(ii).line_x)); % DEM Line
%     for jj = 1:length(T(ii).line_x)
%         dist = sqrt((T(ii).line_x(jj)-t_x).^2 + (T(ii).line_y(jj)-t_y).^2);
%         [~, I] = min(dist);
%         T(ii).line_z(jj) = t_z(I);
%     end
%     clear I jj
%
%     % Reverse transect if needed, want ocean to be first
%     if T(ii).line_z(1) > T(ii).line_z(end)
%         T(ii).line_x = fliplr(T(ii).line_x);
%         T(ii).line_y = fliplr(T(ii).line_y);
%         T(ii).line_z = fliplr(T(ii).line_z);
%     end
%     %     % Smooth the elevation data along transect
%     %     T(ii).z_s = smoothn(T(ii).line_z,3);
%
%     % Make Along Transect Grid
%     T(ii).s = sqrt(T(ii).line_x.^2+T(ii).line_y.^2);
%     T(ii).s = T(ii).s - min(T(ii).s);
%
%     clear t_x t_y t_z
%
% end
