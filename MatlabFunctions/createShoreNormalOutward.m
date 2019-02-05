function [xline,yline] = createShoreNormalOutward(m,x_point,y_point,xgrid,ygrid,zgrid,tran_length)
%Calculate orthogonal transect that points outwards from shore

% Inputs:
%   m: slope of normal line
    % x_point: utm of current point in easting
    % y_point: utm of current point in northing 
    % heading: heading of line to calculate outwards direction
    % xgrid, ygrid, zgrid: grids of elevation data in UTM
    % tran_length: length of the final transect that you want
    
    x1 = linspace(x_point,x_point-100,2);
    y1 = m*(x1 - x_point)+y_point;
    step = 10;
    [ line_x1, line_y1 ] = createTransect( x1(1), y1(1), x1(2), y1(2), step ); % Make transect along line
    x2 = linspace(x_point,x_point+100,2);
    y2 = m*(x2 - x_point)+y_point;
    [ line_x2, line_y2 ] = createTransect( x2(1), y2(1), x2(2), y2(2), step );
    
    % Convert normal lines to Lat Lon 
    [ly1,lx1] = utm2deg(line_x1,line_y1,repmat('10 T',length(line_x1),1));
    [ly2,lx2] = utm2deg(line_x2,line_y2,repmat('10 T',length(line_x2),1));
    
    % Set the Line X and Y to the correct heading 
    if I == 2 % line_x1 is pointing offshore
        line_x = line_x2;
        line_y = line_y2;
    elseif I == 1 % line_x2 is pointing offshore
        line_x = line_x1; 
        line_y = line_y1;
    end
    
    
    % Find closest elevation points to the ends of each line segment 
    nanInds = isnan(zgrid);
    tx = xgrid(~nanInds);
    ty = ygrid(~nanInds);
    tz = zgrid(~nanInds);
    % --------- First Segment 
    d1 = zeros(length(line_x1),1);
    for jj = 1:length(line_x1)
        dist = sqrt((line_x1(jj) - tx).^2 + (line_y1(jj) - ty).^2);
        [~,I] = min(dist);
        tempx = tx(I);
        tempy = ty(I);
        [r,c] = find(xgrid == tempx & ygrid == tempy);
        d1(jj) = zgrid(r,c);
    end
    % --------- Second Segment
    d2 = zeros(length(line_x2),1);
    for jj= 1:length(line_x2)
        dist = sqrt((line_x2(jj) - tx).^2 + (line_y2(jj) - ty).^2);
        [~,I] = min(dist);
        tempx = tx(I);
        tempy = ty(I);
        [r,c] = find(xgrid == tempx & ygrid == tempy);
        d2(jj) = zgrid(r,c);
    end
    
    % Find which way is increasing towards land
    [~,I] = max([max(d1),max(d2)]); % find which line segment has the highest elevatio in it
    if I == 2 % line_x1 is pointing offshore
        line_x = line_x1;
        line_y = line_y1;
    elseif I == 1 % line_x2 is pointing offshore
        line_x = line_x2; 
        line_y = line_y2;
    end

    %  Create transect of preferred distance sepcified above by len
    step = 1;
    [ line_x, line_y ] = createTransect( line_x(1), line_y(1), line_x(end), line_y(end), step );
    % Set x,y to be length of desired transect
    xline = line_x(1:tran_length);
    yline = line_y(1:tran_length);
end

