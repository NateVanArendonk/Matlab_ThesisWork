%% Load in variables of interest
addpath C:\Functions_Matlab
clear all
close all
clc

%%
path = 'E:\Abbas\WCRP\Tacoma\Tier3\XBeach\OB_RW\10yrScenario\SLR280\'; % Path to the model runs
nc_filename = '\xboutput.nc'; % Name of output file

dirNames = strcat(path);
xbr = dir(dirNames); % get list of all runs
xbr(1:2) = [];
xbr(end) = [];
% Sort based on transect order
for ii = 1:length(xbr)
    xbr(ii).number = str2num(erase(xbr(ii).name,'Line'));
end
xbr = nestedSortStruct(xbr,{'number'});

for ii = 1:length(xbr)
    if xbr(ii).isdir == 1
        X(ii).name = xbr(ii).name;
        spaceID = strfind(xbr(ii).name,' ');
        X(ii).GridName = X(ii).name;
        X(ii).GridName(spaceID) = [];
        xgrid = sprintf('/%s_XB_x.grd',X(ii).GridName);
        ygrid = sprintf('/%s_XB_y.grd',X(ii).GridName);
        zgrid = sprintf('/%s_XB_z.grd',X(ii).GridName);
        
        % Load Grid Info
        X(ii).x = load([path X(ii).name xgrid]);
        X(ii).y = load([path X(ii).name ygrid]);
        X(ii).z = load([path X(ii).name zgrid]);
        
        % Load wave info
        X(ii).point_zs = squeeze(ncread([path xbr(ii).name nc_filename],'point_zs')); % Runup
        X(ii).point_xz = squeeze(ncread([path xbr(ii).name nc_filename],'point_xz'));
        X(ii).point_yz = squeeze(ncread([path xbr(ii).name nc_filename],'point_yz'));
        X(ii).zs = squeeze(ncread([path xbr(ii).name nc_filename],'zs')); % Water Level
        
        X(ii).Hs = 4*sqrt(var(X(ii).zs(:,1:end-1)')); % Gets rid of anomoly at end
        
        % Make cross shore 'S' transect
        X(ii).s = sqrt(X(ii).x.^2+X(ii).y.^2);
        X(ii).s = X(ii).s - min(X(ii).s);
        
        % Make Cross shore 'S' transect for point Runup gauge data - BROKEN
        X(ii).point_s = sqrt(X(ii).point_xz.^2+X(ii).point_yz.^2);
        X(ii).point_s = X(ii).point_s - min(X(ii).point_s);
    end
end
return

%% Calculate Runup - Stockdon methodology

% Calculate Set up
time = 1:1:length(X(1).point_zs);

% Find Peaks in Run up and plot
for ii = 1:length(X)
    [pks, locs] = findpeaks(X(ii).point_zs,time,'MinPeakDistance',12,'MinPeakWidth',1);
    cdf = sort(pks,'ascend');
    cdfy = linspace(0,1,length(cdf));
    ind = findnearest(0.98,cdfy);
    % Calculate R2
    X(ii).R2 = cdf(ind);
    if isempty(X(ii).R2)
        X(ii).R2 = 0;
    end
    
    X(ii).s_up = nanmean(X(ii).zs,2);
    s_thresh = 0.05;
    diff_vec = X(ii).s_up - X(ii).z';
    sup_ind = find(diff_vec <= s_thresh);
    X(ii).sup_ind = sup_ind(1);
    X(ii).maxSUP = X(ii).s_up(sup_ind(1));
end



%% Find R2 elevation on profile
RX = zeros(length(X),1);
RY = RX;
RS = RX;
RZ = RX;

metric = 'R2';
for ii = 1:length(X)
    switch metric 
        case 'R2'
            
            eleThresh = X(ii).R2(1);
            inds = find(X(ii).z <= eleThresh);
            dif = sqrt((X(ii).z(inds) - eleThresh).^2); % Calculate distnace between elevation and every point on elevation profile 
            [dif,I] = sort(dif,'ascend'); % Find first instance and closest to begining elevation threshold 
            I = I(1);
     
            RX(ii) = X(ii).x(I);
            RY(ii) = X(ii).y(I);
            RS(ii) = X(ii).s(I);
            RZ(ii) = X(ii).R2(1);
        case 'SUP'
            RX(ii) = X(ii).x(X(ii).sup_ind);
            RY(ii) = X(ii).y(X(ii).sup_ind);
            RS(ii) = X(ii).s(X(ii).sup_ind);
            RZ(ii) = X(ii).s_up(X(ii).sup_ind);
    end         
end

%% Load in DEM and subsample to be just owen beach and ruston way

% First load in the DEM
E = load('E:\Abbas\Modeling Resources\PS_DEM\Ruston_Way\RustonWayCONED_DEM.mat');

K = kml2struct(['E:\Abbas\WCRP\Tacoma\Tier3\KML_IN\' 'OB_RW_Mask.kml']);
[K.x,K.y] = deg2utm(K.Lat,K.Lon);
tic
IN = inpolygon(E.x,E.y,K.x,K.y);
toc
E.x = E.x(IN);
E.y = E.y(IN);
E.z = E.z(IN);

% Load in Zone Mask KMLs
mask_fol = 'E:\Abbas\WCRP\Tacoma\Tier3\XBeach\MaskZones\';
masks = dir('E:\Abbas\WCRP\Tacoma\Tier3\XBeach\MaskZones\*.kml');
for mm = 1:length(masks)
    temp = kml2struct([mask_fol masks(mm).name]);
    Z(mm).name = temp.Name;
    Z(mm).lon = temp.Lon;
    Z(mm).lat = temp.Lat;
    [Z(mm).x,Z(mm).y] = deg2utm(temp.Lat,temp.Lon);
    tic
    IN = inpolygon(E.x,E.y,Z(mm).x,Z(mm).y);
    % Subset DEM
    subX = E.x(IN);
    subY = E.y(IN);
    subZ = E.z(IN);
    Z(mm).mx = min(subX):1:max(subX);
    Z(mm).my = min(subY):1:max(subY);
    [Z(mm).MX,Z(mm).MY] = meshgrid(Z(mm).mx,Z(mm).my);
    Z(mm).MZ = griddata(subX,subY,subZ,Z(mm).MX,Z(mm).MY);
    
    toc
end
clear subX subY subZ E
%% Contour Between Points
% 14 minutes to run
CX = cell(length(X),1);
CY = CX;
a = [];
b = [];
tic

% inds2smooth = [123:130,150:156];
inds2smooth = [];
zone1 = 1:84;
zone2 = 85:134;
zone3 = 135:length(X);
for ii = 2:90%length(X)
    % ---------------------------------------------------------------------
    % R2 values at each transect 
    rz1 = RZ(ii);
    rz2 = RZ(ii-1);
    
    % First Transect XY coordiantes of R2
    x1 = X(ii).x;
    y1 = X(ii).y;
    % Second Tramsect XY coordinates of R2
    x2 = X(ii-1).x;
    y2 = X(ii-1).y;
    
    % Decide which portion DEM to use
    if ismember(ii,zone1)
        mx = Z(1).mx;
        my = Z(1).my;
        MX = Z(1).MX;
        MY = Z(1).MY;
        MZ = Z(1).MZ;
    elseif ismember(ii,zone2)
        mx = Z(2).mx;
        my = Z(2).my;
        MX = Z(2).MX;
        MY = Z(2).MY;
        MZ = Z(2).MZ;
    else
        mx = Z(3).mx;
        my = Z(3).my;
        MX = Z(3).MX;
        MY = Z(3).MY;
        MZ = Z(3).MZ;
    end
    
    % This is used for logic later on, how big to draw a polygon which is
    % used for limiting the contour. 
    % The size of this buffer can slow down this code heavily so try and
    % get away with smallest but biggest polyshape necessary
    % Make sure you have those number of indices available
    % ex. at 25 along ruston, there are ~1800 indices in variable a
    if ii < 25
        ainds = 1:length(a);
    elseif ii > 25
        ainds = length(a)-500:length(a);
    end
    
    % Smooth DEM using Convolution filter for nusiance area 
    if ismember(ii,inds2smooth)
        kernel = ones(2)/2;
        MZ = conv2(MZ,kernel,'same');
    end
    
    % These are the current Transect Lines used for finding lines of
    % intersection
    L1 = [x1;y1];
    L2 = [x2;y2];
    
    % --------------------- First R2 Value Contour ------------------------
    % Go through and find all of the contours created and then which ones
    % go through both of our transects
    P1 = []; % First Transect and First R2
    P3 = []; % Second Transect and First R2
    cc = contourc(mx,my,MZ,[rz1 rz1]); % Make first contour
    cc(:,1) = []; % Get rid of beginning cells that are meaningless for this 
    cx1 = cc(1,:); % Grab X's
    cy1 = cc(2,:); % Grab Y's
    bID = find(cx1 < 1000);% Find breaks in the contour
    if ~isempty(bID) % if there are breaks find the longest continuous contour between transect
        % Loop through and grab indice brackets for contour groups
        inds = zeros(length(bID)+1,2);
        dvec = zeros(length(bID)+1,1);
        for nn = 1:length(bID)+1 % Populate with contour groups 
            if nn == 1
                inds(1,1) = 1;
                inds(1,2) = bID(nn)-1;
                dvec(1) = inds(1,2)-inds(1,1);
            elseif nn == length(bID)+1
                inds(nn,1) = bID(nn-1)+1;
                inds(nn,2) = length(cx1);
                dvec(nn) = inds(nn,2)-inds(nn,1);
            else
                inds(nn,1) = bID(nn-1)+1;
                inds(nn,2) = bID(nn)-1;
                dvec(nn) = inds(nn,2)-inds(nn,1);
            end
        end
        [~,II] = sort(dvec,'descend'); % sort the contour groups with longest continuous contours first 
        inds = inds(II,:); %Grab the top 5 and ditch the rest 
        if length(inds) > 5
            temp_inds = inds(1:5,:);
            clear inds;
            inds = temp_inds;
        end
        % loop through and find the scenario where both ends of contour crosses both transects
        ind2use = 1;
        while isempty(P1) || isempty(P3) 
            block = inds(ind2use,1):inds(ind2use,2);
            tx1 = cx1(block);
            ty1 = cy1(block);
            
            C1 = [tx1;ty1];
            P1 = InterX(L1,C1); % First Transect and First Contour End
            P3 = InterX(L2,C1); % Second Trcansect and First Contour Beginning
            % If these variables are not empty, then the contour corsses
            % both points 
            ind2use = ind2use+1;
        end
    else % Otherwise if it is a single contour, grob the intersection of it with the transects 
        tx1 = cx1;
        ty1 = cy1;
        C1 = [tx1;ty1];
        P1 = InterX(L1,C1); % First Transect and First Contour End
        %----------------
        P3 = InterX(L2,C1); % Second Transect and First Contour Beginning
    end
    
    % --------------------- Second R2 Value Contour -----------------------
    P2 = []; % Second R2 contour with first transect
    P4 = []; % Second R2 contour with second transect
    cc = contourc(mx,my,MZ,[rz2 rz2]); % Make contour
    cc(:,1) = []; % Get rid of first part of data which is meaningless
    cx2 = cc(1,:); % Seperate in to X and Y
    cy2 = cc(2,:); 
    bID = find(cx2 < 1000);% Find breaks in the contour
    if ~isempty(bID)
        % Loop through and grab indice brackets for contour groups
        inds = zeros(length(bID)+1,2);
        dvec = zeros(length(bID)+1,1);
        for nn = 1:length(bID)+1 % Populate indice groups 
            if nn == 1
                inds(1,1) = 1;
                inds(1,2) = bID(nn)-1;
                dvec(1) = inds(1,2)-inds(1,1);
            elseif nn == length(bID)+1
                inds(nn,1) = bID(nn-1)+1;
                inds(nn,2) = length(cx2);
                dvec(nn) = inds(nn,2)-inds(nn,1);
            else
                inds(nn,1) = bID(nn-1)+1;
                inds(nn,2) = bID(nn)-1;
                dvec(nn) = inds(nn,2)-inds(nn,1);
            end
        end
        [~,II] = sort(dvec,'descend'); % Sort by largest continuous contour groups 
        inds = inds(II,:);
        if length(inds) > 5 % Grab the top 5 contours and ditch the rest 
            temp_inds = inds(1:5,:);
            clear inds;
            inds = temp_inds;
        end
        ind2use = 1;
        while isempty(P2) || isempty(P4) % While these are empty keep trying the contours 
            block = inds(ind2use,1):inds(ind2use,2);
            tx2 = cx2(block);
            ty2 = cy2(block);
            
            C2 = [tx2;ty2];
            P2 = InterX(L1,C2); % First Transect and Second Contour End
            P4 = InterX(L2,C2); % Second Transect and Second Contour Beginning
            ind2use = ind2use+1;
        end
    else % otherwise if there is one contour just use taht 
        tx2 = cx2;
        ty2 = cy2;
        C2 = [tx2;ty2];
        P2 = InterX(L1,C2); % First Transect and First Contour End
        P4 = InterX(L2,C2); % Second Transect and First Contour Beginning
    end
    
    
    % -------------------------------------------------------------------------
    % Heavy logic below - very confusing
    % I1 and I2 have same logic an they are assumed to be limit of new
    % build of alongshore contour
    % I3 and I4 have same logic and they are assumed to be start of new
    % build of alongshore contour 
    % Good Luck 
    [rr,cc] = size(P1); %If there are multiple crossing of the contour on the transect, go through the logic 
    if cc > 1
        px = P1(1,:); % grab x and y variables of crossings
        py = P1(2,:);
        distances = ones(length(px),1); % Loop through and see if the point is within the alongshore contour already, a bit unnecessary for this point, but a check 
        for zz = 1:length(px)
            dist = sqrt((px(zz) - a).^2 + (py(zz) - b).^2); 
            dmin = min(dist);
            if dmin <= 1 % if the point is within 5m of the current alongshore contour don't draw that point 
                distances(zz) = 0;
            end
        end
        distances = logical(distances);
        px = px(distances); % get rid of points that are close
        py = py(distances);
        if isempty(px) % if all points are close 
            % Find which segemnt is the least in the currently drawn
            % polygon and then draw that segment.
            px = P1(1,:);
            py = P1(2,:);
            tot = zeros(length(px),1);
            for ll = 1:length(px)
                dist = sqrt((px(ll)-tx1).^2 + (py(ll)-ty1).^2);
                [~,nd] = min(dist);
                % Make a polygon around the current alongshore transect
                utmzone = repmat('10 T',length(ainds),1);
                [lat,lon] = utm2deg(a(ainds),b(ainds),utmzone);
                bufwidth = 0.00001;
                [latb,lonb] = bufferm(lat,lon,bufwidth);
                [xb,yb] = deg2utm(latb,lonb);
                inds = nd-100:nd+100;
                inp = inpolygon(tx1(inds),ty1(inds),xb,yb);
                tot(ll) = sum(inp);
            end
            [~,imin] = min(tot);
            px = px(imin);
            py = py(imin);
        end
        contourInd = zeros(length(px),1); % loop through and grab the indice on the contour that intersects the point
        flag = ones(length(px),1);
        for pp = 1:length(px)
            dist = sqrt((px(pp)-tx1).^2 + (py(pp)-ty1).^2);
            [~,contourInd(pp)] = min(dist);
            temp_inds = contourInd(pp)-1000:contourInd(pp)+1000;
            if temp_inds(end) > length(tx1)
                temp_inds = contourInd(pp)-1000:length(tx1);
            end
            if temp_inds(1) < 1
                temp_inds = 1:contourInd(pp)+1000;
            end
            % Make a a 'box' of the current transects and see which direction is mostly in the polygon box
            indsIn = inpolygon(tx1(temp_inds),ty1(temp_inds),[x2(1),x2(end),x1(end),x1(1), x2(1)],[y2(1),y2(end),y1(end),y1(1), y2(1)]); 
            if sum(indsIn) > 5 % Give a 5 point buffer
                flag(pp) = 0;
            end
        end
        flag = logical(flag);
        if sum(flag) == 0% If they all fail the logic, then check and see which one is closest to calculated R2 value 
            px = P1(1,:);
            py = P1(2,:);
            dist = sqrt((px-RX(ii)).^2 + (py-RY(ii)).^2); % distance to R2 xy point
            [~,imin] = min(dist);
            px = px(imin);
            py = py(imin);
            dist = sqrt((px-tx1).^2 + (py-ty1).^2);
            [~,imin] = min(dist);
            contourInd = imin;
        else
            contourInd = contourInd(flag);
        end
        if numel(contourInd) == 1
            I1 = contourInd;
        else
            disp('Broken Code come fix I1')
            break
        end
    else % otherwise just find the intersection if there is only one contour
        dist = sqrt((P1(1)-tx1).^2 + (P1(2)-ty1).^2);
        [~, I1] = min(dist); % keep this min bcs it's assuming only 1 crossing
    end
    
    
    
    
    
    
    [rr,cc] = size(P2); %If there are multiple crossing do some fun logic
    if cc > 1
        px = P2(1,:); % grab x and y variables of crossings
        py = P2(2,:);
        distances = ones(length(px),1);
        for zz = 1:length(px)
            dist = sqrt((px(zz) - a(end)).^2 + (py(zz) - b(end)).^2); % This is checking to see if the starting point is already on the alongshore contour
            dmin = min(dist);
            if dmin <= 1 % if the point is within 5m of the current alongshore contour, don't draw it
                distances(zz) = 0;
            end
        end
        distances = logical(distances);
        px = px(distances); % get rid of points that are close
        py = py(distances);
        if isempty(px) % if both are close to it then do next bit of logic
            % Find which segemnt is the least in the currently drawn
            % polygon and then draw that segment.
            px = P2(1,:);
            py = P2(2,:);
            tot = zeros(length(px),1);
            for ll = 1:length(px)
                dist = sqrt((px(ll)-tx2).^2 + (py(ll)-ty2).^2);
                [~,nd] = min(dist);
                % Make a polygon around the current alongshore transect
                utmzone = repmat('10 T',length(a(ainds)),1);
                [lat,lon] = utm2deg(a(ainds),b(ainds),utmzone);
                bufwidth = 0.00001;
                [latb,lonb] = bufferm(lat,lon,bufwidth);
                [xb,yb] = deg2utm(latb,lonb);
                inds = nd-100:nd+100;
                inp = inpolygon(tx2(inds),ty2(inds),xb,yb);
                tot(ll) = sum(inp);
            end
            [~,imin] = min(tot);
            px = px(imin);
            py = py(imin);
        end
        contourInd = zeros(length(px),1); % loop through and grab the indice on the contour that intersects the point
        flag = ones(length(px),1);
        for pp = 1:length(px)
            dist = sqrt((px(pp)-tx2).^2 + (py(pp)-ty2).^2);
            [~,contourInd(pp)] = min(dist);
            temp_inds = contourInd(pp)-1000:contourInd(pp)+1000;
            if temp_inds(end) > length(tx2)
                temp_inds = contourInd(pp)-1000:length(tx2);
            end
            if temp_inds(1) < 1
                temp_inds = 1:contourInd(pp)+1000;
            end
            indsIn = inpolygon(tx2(temp_inds),ty2(temp_inds),[x2(1),x2(end),x1(end),x1(1)],[y2(1),y2(end),y1(end),y1(1)]);
            if sum(indsIn) > 5 % Give a 5 point buffer
                flag(pp) = 0;
            end
        end
        flag = logical(flag);
        if sum(flag) == 0% If they all fail the logic, then check and see which one is closest to calculated R2 value 
            px = P2(1,:);
            py = P2(2,:);
            dist = sqrt((px-RX(ii)).^2 + (py-RY(ii)).^2); % distance to R2 xy point
            [~,imin] = min(dist);
            px = px(imin);
            py = py(imin);
            dist = sqrt((px-tx2).^2 + (py-ty2).^2);
            [~,imin] = min(dist);
            contourInd = imin;
        else
            contourInd = contourInd(flag);
        end
        if numel(contourInd) == 1
            I2 = contourInd;
        else
            disp('Broken Code come fix I2')
            break
        end
    else % otherwise just find the intersection if there is only one contour
        dist = sqrt((P2(1)-tx2).^2 + (P2(2)-ty2).^2);
        [~, I2] = min(dist); % keep this min bcs it's assuming only 1 crossing
    end
    
    
     
    [rr,cc] = size(P3); %If there are multiple crossing do some fun logic
    if cc > 1
        px = P3(1,:); % grab x and y variables of crossings
        py = P3(2,:);
        distances = ones(length(px),1);
        for zz = 1:length(px)
            dist = sqrt((px(zz) - a(end)).^2 + (py(zz) - b(end)).^2); % This is checking to see if the starting point is already on the alongshore contour
            dmin = min(dist);
            if dmin <= 5 % if the point is within 5m of the current alongshore contour, don't draw it
                distances(zz) = 1;
            else
                distances(zz) = 0;
            end
        end
        distances = logical(distances);
        px = px(distances); % get rid of points that are close
        py = py(distances);
        if isempty(px) % if both are close to it then do next bit of logic
            % Find which segemnt is the least in the currently drawn
            % polygon and then draw that segment.
            
            px = P3(1,:);
            py = P3(2,:);
            
            dist = sqrt((px-RX(ii-1)).^2 + (py-RY(ii-1)).^2); % distance to R2 xy point
            [~,imin] = min(dist);
            px = px(imin);
            py = py(imin);
            dist = sqrt((px-tx1).^2 + (py-ty1).^2);
            [~,imin] = min(dist);
            contourInd = imin;
            
%             % Make a polygon around the current alongshore transect
%             utmzone = repmat('10 T',length(a(ainds)),1);
%             [lat,lon] = utm2deg(a(ainds),b(ainds),utmzone);
%             bufwidth = 0.00001;
%             [latb,lonb] = bufferm(lat,lon,bufwidth);
%             [xb,yb] = deg2utm(latb,lonb);
%             px = P3(1,:);
%             py = P3(2,:);
%             tot = zeros(length(px),1);
%             for ll = 1:length(px)
%                 dist = sqrt((px(ll)-tx1).^2 + (py(ll)-ty1).^2); % find point on contour of ith intersection
%                 [~,nd] = min(dist);
%                 inds = nd-100:nd+100;
%                 inp = inpolygon(tx1(inds),ty1(inds),xb,yb);
%                 tot(ll) = sum(inp);
%             end
%             [~,imin] = min(tot);
%             px = px(imin);
%             py = py(imin);
%             dist = sqrt((px-tx1).^2 + (py-ty1).^2);
%             [~,contourInd] = min(dist);
        elseif numel(px) > 1
            contourInd = zeros(length(px),1); % loop through and grab the indice on the contour that intersects the point
            %         flag = ones(length(px),1);
            utmzone = repmat('10 T',length(a(ainds)),1);
            [lat,lon] = utm2deg(a(ainds),b(ainds),utmzone);
            bufwidth = 0.000012;
            [latb,lonb] = bufferm(lat,lon,bufwidth);
            [xb,yb] = deg2utm(latb,lonb);
            tot = zeros(length(px),1);
            for pp = 1:length(px)
                dist = sqrt((px(pp)-tx1).^2 + (py(pp)-ty1).^2);
                [~,contourInd(pp)] = min(dist);
                cstart = contourInd(pp)-500;
                cend = contourInd(pp) + 500;
                if cstart < 1
                    cstart = 1;
                end
                if cend > length(tx1)
                    cend = length(tx1);
                end
                temp_inds = cstart:cend;
                indsIn = inpolygon(tx1(temp_inds),ty1(temp_inds),xb,yb);
                tot(pp) = sum(indsIn);
            end
            [~,imax] = max(tot);
            contourInd = contourInd(imax);
        else
            dist = sqrt((px-tx1).^2 + (py-ty1).^2);
            [~,contourInd] = min(dist);
        end
        if numel(contourInd) == 1
            I3 = contourInd;
        else
            disp('Broken Code come fix I3')
            break
        end
    else % otherwise just find the intersection if there is only one contour
        dist = sqrt((P3(1)-tx1).^2 + (P3(2)-ty1).^2);
        [~, I3] = min(dist); % keep this min bcs it's assuming only 1 crossing
    end
    
    
    
    
    
    
    [rr,cc] = size(P4); %If there are multiple crossing do some fun logic
    if cc > 1
        px = P4(1,:); % grab x and y variables of crossings
        py = P4(2,:);
        distances = ones(length(px),1);
        for zz = 1:length(px)
            dist = sqrt((px(zz) - a(end)).^2 + (py(zz) - b(end)).^2); % This is checking to see if the starting point is already on the alongshore contour
            dmin = min(dist);
            if dmin <= 5 % if the point is within 5m of the current alongshore contour, don't draw it
                distances(zz) = 1;
            else
                distances(zz) = 0;
            end
        end
        distances = logical(distances);
        px = px(distances); % get rid of points that are close
        py = py(distances);
        if isempty(px) % if both are close to it then do next bit of logic
            % Find which segemnt is the least in the currently drawn
            % polygon and then draw that segment.
            px = P4(1,:);
            py = P4(2,:);
            
            dist = sqrt((px-RX(ii-1)).^2 + (py-RY(ii-1)).^2); % distance to R2 xy point
            [~,imin] = min(dist);
            px = px(imin);
            py = py(imin);
            dist = sqrt((px-tx2).^2 + (py-ty2).^2);
            [~,imin] = min(dist);
            contourInd = imin;
%             % Make a polygon around the current alongshore transect
%             utmzone = repmat('10 T',length(a(ainds)),1);
%             [lat,lon] = utm2deg(a(ainds),b(ainds),utmzone);
%             bufwidth = 0.00001;
%             [latb,lonb] = bufferm(lat,lon,bufwidth);
%             [xb,yb] = deg2utm(latb,lonb);
%             px = P4(1,:);
%             py = P4(2,:);
%             tot = zeros(length(px),1);
%             for ll = 1:length(px)
%                 dist = sqrt((px(ll)-tx2).^2 + (py(ll)-ty2).^2); % find point on contour of ith intersection
%                 [~,nd] = min(dist);
%                 inds = nd-100:nd+100;
%                 inp = inpolygon(tx2(inds),ty2(inds),xb,yb);
%                 tot(ll) = sum(inp);
%             end
%             [~,imin] = min(tot);
%             px = px(imin);
%             py = py(imin);
%             dist = sqrt((px-tx2).^2 + (py-ty2).^2);
%             [~,contourInd] = min(dist);
%             
        elseif numel(px) > 1
            contourInd = zeros(length(px),1); % loop through and grab the indice on the contour that intersects the point
            %         flag = ones(length(px),1);
            utmzone = repmat('10 T',length(a(ainds)),1);
            [lat,lon] = utm2deg(a(ainds),b(ainds),utmzone);
            bufwidth = 0.000012;
            [latb,lonb] = bufferm(lat,lon,bufwidth);
            [xb,yb] = deg2utm(latb,lonb);
            tot = zeros(length(px),1);
            for pp = 1:length(px)
                dist = sqrt((px(pp)-tx2).^2 + (py(pp)-ty2).^2);
                [~,contourInd(pp)] = min(dist);
                cstart = contourInd(pp)-500;
                cend = contourInd(pp) + 500;
                if cstart < 1
                    cstart = 1;
                end
                if cend > length(tx2)
                    cend = length(tx2);
                end
                temp_inds = cstart:cend;
                indsIn = inpolygon(tx2(temp_inds),ty2(temp_inds),xb,yb);
                tot(pp) = sum(indsIn);
            end
            [~,imax] = max(tot);
            contourInd = contourInd(imax);
        else
            dist = sqrt((px-tx2).^2 + (py-ty2).^2);
            [~,contourInd] = min(dist);
        end
        if numel(contourInd) == 1
            I4 = contourInd;
        else
            disp('Broken Code come fix I4')
            break
        end
    else % otherwise just find the intersection if there is only one contour
        dist = sqrt((P4(1)-tx2).^2 + (P4(2)-ty2).^2);
        [~, I4] = min(dist); % keep this min bcs it's assuming only 1 crossing
    end
    
    
    
    % Subsample contour now - ix1 and iy1 should be the same
    % Just want portin of contour between transects
    
    inds1 = [I1,I3];
    inds2 = [I2,I4];
    inds1 = sort(inds1);
    inds2 = sort(inds2);
    tx1 = tx1(inds1(1):inds1(end));
    ty1 = ty1(inds1(1):inds1(end));
    tx2 = tx2(inds2(1):inds2(end));
    ty2 = ty2(inds2(1):inds2(end));
    
    len1 = length(tx1);
    len2 = length(tx2);
    lens = [len1;len2];
    
    % If the two lengths don't equal eachother, do this little hack to make
    % them - Don't think this method is correct but will do for
    % right now
    if len1 ~= len2
        [~,II] = max(lens);
        if II == 1
            fx = [];
            fy = [];
            for nn = 1:length(tx1)-1
                step = .001; % 100 meters,  dx = res_in_meters/73000, Note I used 91500 its half between 110000 and 73000
                [tempx, tempy ] = createTransect( tx1(nn), ty1(nn), tx1(nn+1), ty1(nn+1), step );
                fx = cat(1,fx,tempx');
                fy = cat(1,fy,tempy');
            end
            inds = 1:round(length(fx)/len2):length(fx);
            tx1 = fx(inds);
            ty1 = fy(inds);
%             if length(tx1) ~= len2
%                 len1 = length(tx1);
%                 len2 = length(tx2);
%                 toCut = abs(len2 - length(tx1));
%                 tx1(1:toCut) = [];
%                 ty1(1:toCut) = [];
%             end
%             tx1 = tx1';
%             ty1 = ty1';
        else
            fx = [];
            fy = [];
            for nn = 1:length(tx2)-1
                step = .001; % 
                [tempx, tempy ] = createTransect( tx2(nn), ty2(nn), tx2(nn+1), ty2(nn+1), step );
                fx = cat(1,fx,tempx');
                fy = cat(1,fy,tempy');
            end
            inds = 1:round(length(fx)/len1):length(fx);
            tx2 = fx(inds);
            ty2 = fy(inds);
%             if length(tx2) ~= len1
%                 % Trim up the longer one 
%                 len1 = length(tx1);
%                 len2 = length(tx2);
%                 t2cut = [len1;len2];
%                 [~,II] = max(t2cut);
%                 % Find how many indices to cut 
%                 toCut = abs(len1 - length(tx2));
%                 tx2(1:toCut) = [];
%                 ty2(1:toCut) = [];
%             end
%             tx2 = tx2';
%             ty2 = ty2';
        end
        if length(tx2) ~= length(tx1) % double check to see if they're the same length, if not, just snip off the difference -_-
            len1 = length(tx1);
            len2 = length(tx2);
            t2cut = [len1;len2];
            [~,II] = max(t2cut);
            % Find how many indices to cut
            toCut = abs(len1 - length(tx2));
            if II == 1
                tx1(1:toCut) = [];
                ty1(1:toCut) = [];
                tx1 = tx1';
                ty1 = ty1';
            else
                tx2(1:toCut) = [];
                ty2(1:toCut) = [];
                tx2 = tx2';
                ty2 = ty2';
            end
        end
                
            
    end
    alpha = 0:1/(length(tx1)-1):1;
    % Make sure dimensions of tx's and ty's are correct
    [rr,~] = size(tx1);
    if rr ~= 1
        tx1 = tx1';
        ty1 = ty1';
    end
    [rr,~] = size(tx2);
    if rr ~= 1
        tx2 = tx2';
        ty2 = ty2';
    end
    
    CX{ii} = tx1.*alpha + tx2.*(1-alpha);
    CY{ii} = ty1.*alpha + ty2.*(1-alpha);
    a = cat(2,a,CX{ii});
    b = cat(2,b,CY{ii});
    
    %                 clf
    %                 pcolor(MX,MY,MZ)
    %                 shading interp
    %                 hold on
    %                 plot(x1,y1,'k')
    %                 plot(x2,y2,'k')
    %                 plot(tx1,ty1)
    %                 plot(tx1,ty1,'k')
    %                 tit = sprintf('Transect %d',ii);
    %                 title(tit)
    %                 pause(0.3)
    
    if rem(ii,10) == 0
        fprintf('Completed %d out of %d\n',ii,length(X))
    end
end
toc








% Make KML 

% a = [];
% b = [];
% for ii = 2:length(CX)
%     a = cat(2,a,CX{ii});
%     b = cat(2,b,CY{ii});
% end


utmzone = repmat('10 T',length(a),1);
[lat,lon] = utm2deg(a,b,utmzone);

% Write KML
% kmlwriteline('OBRW_10yr_280FTSLR.kml',lat,lon);
kmlwriteline('OBRW_SUP_10yr_SLR280.kml',lat,lon);
fprintf('KML Line Saved\n')
