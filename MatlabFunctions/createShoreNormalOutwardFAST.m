function [xline,yline] = createShoreNormalOutwardFAST(m,x_point,y_point,heading,tran_length)
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
    
    % Calculate heading of each line 
    [~,h1] = distance(ly1(1),lx1(1),ly1(end),lx1(end)); 
    [~,h2] = distance(ly2(1),lx2(1),ly2(end),lx2(end));
    
    % Determine which heading is closest to 90 off the one we prescribe
    headingWant = wrapTo360(heading - 90);
    
    % Take the absolute difference between the heading we want and the two
    % options for heading 
    d1 = abs(headingWant - h1);
    d2 = abs(headingWant - h2);
    
    % Find which one has the smallest difference 
    [~,I] = min([d1,d2]);
    
    % Set the Line X and Y to the correct heading 
    if I == 2 % line_x1 is pointing offshore
        line_x = line_x2;
        line_y = line_y2;
    elseif I == 1 % line_x2 is pointing offshore
        line_x = line_x1; 
        line_y = line_y1;
    end
    
    %  Create transect of preferred distance sepcified above by len
    step = 1;
    [ line_x, line_y ] = createTransect( line_x(1), line_y(1), line_x(end), line_y(end), step );
    % Set x,y to be length of desired transect
    xline = line_x(1:tran_length);
    yline = line_y(1:tran_length);
end

