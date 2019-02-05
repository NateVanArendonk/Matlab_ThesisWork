clearvars
clc
clf

% Load in shoreline KML
kml_fol = 'C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\KML\PugetSoundShorline\';
K = kml2struct([kml_fol 'RustonWay.kml']);

% Load in Bathy/DEM
D = load('E:\Abbas\Modeling Resources\PS_DEM\Ruston_Way\RustonWayCONED_DEM.mat');

% ----------------- USER INPUT ------------------------------------
spacing = 100; %[m] Spacing between transects
len = 300; %[m] Length of transect
% -----------------------------------------------------------------

% Convert to KML to  UTM 
[K.x,K.y] = deg2utm(K.Lat,K.Lon);
line(K.x,K.y)
hold on

% Make polygon to confirm what line to use 
% bay_poly = polyshape([K.x' K.x(end) K.x(1)],[K.y' 5.392*10^6 5.392*10^6]);
% plot(bay_poly)

% Create finely spaced alongshore transect to be derefined later
fx_e = [];
fy_e = [];
for nn = 1:length(K.x)-1
    step = .5;
    [ line_x, line_y ] = createTransect( K.x(nn), K.y(nn), K.x(nn+1), K.y(nn+1), step );
    fx_e = cat(1,fx_e,line_x');
    fy_e = cat(1,fy_e,line_y');
end
% Get rid of any nans
nan_inds = isnan(fx_e);
fx_e(isnan(fx_e)) = [];
fy_e(isnan(fy_e)) = [];
K.x = fx_e;
K.y = fy_e;
clf
line(K.x,K.y)
hold on 

inds = 0:spacing:length(K.x);
inds(1) = 2;

X = zeros(length(inds),2);


% Loop through and calculate perpendicular line
for ii = 1:length(inds)
    % Get current point, and previous point for calculation of slope of normal line
    x_p = [K.x(inds(ii)-1), K.x(inds(ii))];
    y_p = [K.y(inds(ii)-1), K.y(inds(ii))];
    
    % Create finely spaced line between the two points
    x_e = [];
    y_e = x_e;
    step = 0.01;
    [ line_x, line_y ] = createTransect( x_p(1), y_p(1), x_p(2), y_p(2), step );
    x_e = cat(1,x_e,line_x');
    y_e = cat(1,y_e,line_y');
    its = round(length(x_e)-1);
    if its >= 5
        % Loop through and grab average slope between first point and
        % subsequent point in line between transects
        slope = zeros(its,1);
        for jj = 2:its
            slope(jj,1) = (y_e(jj) - y_e(jj-1))/(x_e(jj) - x_e(jj-1));
        end
        % Get rid of any nans or zeros
        slope(isnan(slope)) = [];
        slope(slope == 0) = [];
        % Get average slope of all points 
        slope = mean(slope);
        % Calculate slope of normal line
        m = -1/slope;
        
        % Create Orthogonal Line extending beyond the point of interest in
        % both directions 
        x1 = linspace(x_p(1),x_p(1)-100,2);
        y1 = m*(x1 - x_p(1))+y_p(1);
        step = 10;
        [ line_x1, line_y1 ] = createTransect( x1(1), y1(1), x1(2), y1(2), step ); % Make transect along line
        x2 = linspace(x_p(1),x_p(1)+100,2);
        y2 = m*(x2 - x_p(1))+y_p(1);
        [ line_x2, line_y2 ] = createTransect( x2(1), y2(1), x2(2), y2(2), step );
        
        % Grab bathy along each part of the extended line
        % Extract subset of depths that covers transect
        t1_inds = D.x >= min(line_x1) & D.x <= max(line_x1) & ...
            D.y >= min(line_y1) & D.y <= max(line_y1);
        t1_x = D.x(t1_inds);
        t1_y = D.y(t1_inds);
        t1_z = D.z(t1_inds);
        t2_inds = D.x >= min(line_x2) & D.x <= max(line_x2) & ...
            D.y >= min(line_y2) & D.y <= max(line_y2);
        t2_x = D.x(t2_inds);
        t2_y = D.y(t2_inds);
        t2_z = D.z(t2_inds);
        
        
        % Find nearest depth for each location on transect
        line_z1 = zeros(size(line_x1)); % DEM Line
        line_z2 = line_z1;
        for jj = 1:length(line_x1)
            dist = (line_x1(jj)-t1_x).^2 + (line_y1(jj)-t1_y).^2;
            [~, I] = min(dist);
            line_z1(jj) = t1_z(I);
            dist = (line_x2(jj)-t2_x).^2 + (line_y2(jj)-t2_y).^2;
            [~, I] = min(dist);
            line_z2(jj) = t2_z(I);
        end
        
        % Find which way is increasing towards land based on which
        % line_z(end) value is larger and make a super long transect that
        % way
        if line_z1(end) > line_z2(end)      
            x = linspace(x_p(1),x_p(1)-500,2);
            y = m*(x - x_p(1))+y_p(1);       
        else
            x = linspace(x_p(1),x_p(1)+500,2);
            y = m*(x - x_p(1))+y_p(1);
        end
        
%         Create transect of preferred distance sepcified above by len
        tx_e = [];
        ty_e = tx_e;
        step = 1;
        [ line_x, line_y ] = createTransect( x(1), y(1), x(2), y(2), step );
        tx_e = cat(1,tx_e,line_x');
        ty_e = cat(1,ty_e,line_y');
        % Set x,y to be length of desired transect
        x = tx_e(1:len);
        y = ty_e(1:len);
        X(ii,:) = [x(1),x(end)];
        Y(ii,:) = [y(1),y(end)];
        plot(x,y,'k')
    end
    
    plotting = 0; % Only turn on for debugging purposes - slow it down
    if plotting
        hold on
        plot(x,y,'k')
        axis equal
    end
end
XU = X;
YU = Y; 
% Convert back to lat lon for making KML
zeroInds = find(X == 0);
utmzone = repmat('10 U',length(X),1);
[Y(:,1),X(:,1)] = utm2deg(X(:,1),Y(:,1),utmzone);
[Y(:,2),X(:,2)] = utm2deg(X(:,2),Y(:,2),utmzone);

% Sub sample to be desired distance between transects  
% X = X(1:spacing*2:end,:);
% Y = Y(1:spacing*2:end,:);

save('Ruston_Orthog_Transect','X','Y','XU','YU')

%% Write KML file
fid = fopen('Ruston_Orthog_Transects.kml','w');
fprintf(fid,'<?xml version="1.0" encoding="utf-8"?>\n');
fprintf(fid,'<kml xmlns="http://www.opengis.net/kml/2.2">\n');
fprintf(fid,'   <Document>\n');
fprintf(fid,'      <name>Transect_lines</name>\n');
for ii = 1:length(X)
    if Y(ii,1) ~= 0 
        fprintf(fid,'      <Placemark>\n');
        fprintf(fid,'         <Snippet maxLines="0"> </Snippet>\n');
        fprintf(fid,'         <description> </description>\n');
        fprintf(fid,'         <name>Line %d</name>\n',ii);
        fprintf(fid,'         <LineString>\n');
        fprintf(fid,'            <coordinates> %f,%f,0 %f,%f,0</coordinates>\n',X(ii,1),Y(ii,1),X(ii,2),Y(ii,2));
        fprintf(fid,'         </LineString>\n');
        fprintf(fid,'      </Placemark>\n');
    end
end
fprintf(fid,   '</Document>\n');
fprintf(fid,'</kml>\n')
fclose(fid);
fclose('all')
%%
fid = fopen('OB_RW_Orthog_Transects.kml','w');
fprintf(fid,'<?xml version="1.0" encoding="utf-8"?>\n');
fprintf(fid,'<kml xmlns="http://www.opengis.net/kml/2.2">\n');
fprintf(fid,'   <Document>\n');
fprintf(fid,'      <name>Transect_lines</name>\n');
for ii = 1:length(T)
    if T(ii).lat(1) ~= 0 
        fprintf(fid,'      <Placemark>\n');
        fprintf(fid,'         <Snippet maxLines="0"> </Snippet>\n');
        fprintf(fid,'         <description> </description>\n');
        fprintf(fid,'         <name>Line %d</name>\n',ii);
        fprintf(fid,'         <LineString>\n');
        fprintf(fid,'            <coordinates> %f,%f,0 %f,%f,0</coordinates>\n',T(ii).lon(1),T(ii).lat(1),T(ii).lon(2),T(ii).lat(2));
        fprintf(fid,'         </LineString>\n');
        fprintf(fid,'      </Placemark>\n');
    end
end
fprintf(fid,   '</Document>\n');
fprintf(fid,'</kml>\n')
fclose(fid);
fclose('all')