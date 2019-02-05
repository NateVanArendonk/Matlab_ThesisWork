function outAngle = CalculateWaveAttackAngle(X,Y)
% Given Lat and Lon coordinates of transect, calculate angle of wave attack
% or heading or 'Wave Coming From' to be used for XBeach
%   X = x coordinates of transect - Lon
%   Y = y coordinates of transect - Lat
%   Z = elevation coordinates along transect going X and Y - Unnecessary

% if Z(1) > Z(end) % Want First point to indicate land for correct 'heading' calculation 
%     lx1 = X(1);
%     ly1 = Y(1);
%     lx2 = X(end);
%     ly2 = Y(end);
% else
%     lx1 = X(end);
%     ly1 = Y(end);
%     lx2 = X(1);
%     ly2 = Y(1);
% end


% % Some stupid fucking error if this doesn't happen
% restoredefaultpath
% rehash toolboxcache
outAngle = azimuth('gc',Y(1),X(1),Y(2),X(2));
end

