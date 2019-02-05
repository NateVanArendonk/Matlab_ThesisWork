function [U,V] = wind2UV(direct, deg, type, spd)

% Variable direct:
    % direct is in reference to the notation of the wind.  Is the wind
    % being given as heading towards or coming from.  Typically it the wind
    % is given as 'coming from' the south for example - a southerly wind is
    % from the south.  NNRP winds are given as 'From The' notation.  
    % If the wind is from, input: 'from'  otherwise input: 'to'
    
% Varabile deg:
    % This is the wind direction given as an integer 
    
% Variable type:
    % This refers to if the wind is given as cartesian or polar (compass)
    % If the winds are cartesian input: 'cartesian'
    % If the winds are polar input: 'compass'
    
% Variable spd:
    % Wind speed given as an integer

    
% Note - I wrote this script a long time ago and I'm not entirely sure that
% it is working properly but it seems to pass test cases so maybe it does
% do the job 
    

% ex. wind2UV('from', 315, 'cartesian', 10)
type = upper(type);
direct = upper(direct);

if ~length(strfind(type, 'CART')) == 1 % if the wind type is compass
    if strcmpi(direct, 'FROM') % if the winds are given as coming from instead of going
        for ii = 1:length(deg)
            if deg(ii) > 180 && deg(ii) < 360% Convert winds to where they are going
                deg(ii) = deg(ii) - 180;
            elseif deg(ii) == 360
                deg(ii) = 180;
            else
                deg(ii) = deg(ii) + 180;
            end
            % Now convert from compass to cartesian
            deg(ii) = 90 - deg(ii);
        end
    end
else  % Otherwise, if the winds are already in cartesian
    if strcmpi(direct, 'FROM') % if the winds are given as coming from instead of going
        for ii = 1:length(deg)
            if deg(ii) > 180 && deg(ii) < 360
                deg(ii) = deg(ii) - 180;
            elseif deg(ii) == 360
                deg(ii) = 180;
            else
                deg(ii) = deg(ii) + 180;
            end
        end
    end
end
    

% Now break down into U and V components
U = spd.*(cosd(deg));
V = spd.*(sind(deg));




end

    
    
    
    





