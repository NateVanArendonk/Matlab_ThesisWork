%% Load in variables of interest
addpath C:\Functions_Matlab
clear all
close all
clc

% Load in SUP Contour
C = dir('E:\Abbas\WCRP\Tacoma\Tier3\KML_Out\10_yrScenario\SUP\*.kml');
for ii = 1:length(C)
    T(ii) = kml2struct(['E:\Abbas\WCRP\Tacoma\Tier3\KML_Out\10_yrScenario\SUP\' C(ii).name]);
end
for ii = 1:length(T)
    [T(ii).x,T(ii).y] = deg2utm(T(ii).Lat,T(ii).Lon);
    kmlInd = strfind(C(ii).name,'.kml');
    T(ii).Name = C(ii).name(1:kmlInd-1);
end

% Load in All of the Alongshore Transects - to be used later
O = kml2struct('E:\Abbas\WCRP\Tacoma\Tier3\OrthogTran\RW_OB_Transects.kml');
for ii = 1:length(O)
    [O(ii).x,O(ii).y] = deg2utm(O(ii).Lat,O(ii).Lon);
end
D = load('RustonSubset_Gridded_Coned.mat');
%% Load Xbeach Runs
for xx = 1:length(T)
    path = 'E:\Abbas\WCRP\Tacoma\Tier3\XBeach\OB_RW\10yrScenario\'; % Path to the model runs
    
    % Find which run it is and load it 
    ind = strfind(T(xx).Name,'_');
    runName = T(xx).Name(ind(end)+1:length(T(xx).Name));
    path = strcat(path,runName,'\');
    
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
    
    %% Calculate R2 and Sup
    
    time = 1:1:length(X(1).point_zs);
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
        
        % Calculate Setup
        X(ii).s_up = nanmean(X(ii).zs,2);
        s_thresh = 0.05; % 5 cm depth cutoff for setup
        diff_vec = X(ii).s_up - X(ii).z';
        sup_ind = find(diff_vec <= s_thresh);
        X(ii).sup_ind = sup_ind(1);
        X(ii).maxSUP = X(ii).s_up(sup_ind(1));
        
    end
    % Grab only setup vectors and everything
    sx = [];
    sy = [];
    sz = [];
    
    for ii = 1:length(X)
        S(ii).x = X(ii).x(1:X(ii).sup_ind);
        S(ii).y = X(ii).y(1:X(ii).sup_ind);
        S(ii).z = X(ii).z(1:X(ii).sup_ind);
        S(ii).s = X(ii).s(1:X(ii).sup_ind);
        S(ii).sup = X(ii).s_up(1:X(ii).sup_ind);
        sx = horzcat(sx,S(ii).x);
        sy = horzcat(sy,S(ii).y);
        sz = horzcat(sz,S(ii).sup');
    end
    
    %% Load in Gridded Dem
    % First load in the DEM
    % E = load('E:\Abbas\Modeling Resources\PS_DEM\Ruston_Way\RustonWayCONED_DEM.mat');
    % K = kml2struct(['E:\Abbas\WCRP\Tacoma\Tier3\Contour\' 'OB_RW_Mask.kml']);
    % [K.x,K.y] = deg2utm(K.Lat,K.Lon);
    % tic
    % IN = inpolygon(E.x,E.y,K.x,K.y);
    % toc
    % E.x = E.x(IN);
    % E.y = E.y(IN);
    % E.z = E.z(IN);
    %
    % [X,Y,Z] = xyz2grid(E.x,E.y,E.z);

    % Made with XYZ2GRID
    
    
    %% Interp SUP contour on to DEM grid
    SZ = griddata(sx,sy,sz,D.X,D.Y);
    %% Make Contour of LowTide line - Arbitrary Depth
    wl_val = -12; % navd
    cc = contour(D.X,D.Y,D.Z,[wl_val wl_val]);
    cc(:,1) = []; % Get rid of beginning cells that are meaningless for this
    ex = cc(1,:); % Grab X's
    ey = cc(2,:); % Grab Y's
    bID = find(ex < 1000);% Find breaks in the contour
    
    if ~isempty(bID) % if there are breaks find the longest continuous contour between transect
        % Loop through and grab indice brackets for contour groups
        inds2 = zeros(length(bID)+1,2);
        dvec = zeros(length(bID)+1,1);
        for nn = 1:length(bID)+1 % Populate with contour groups
            if nn == 1
                inds2(1,1) = 1;
                inds2(1,2) = bID(nn)-1;
                dvec(1) = inds2(1,2)-inds2(1,1);
            elseif nn == length(bID)+1
                inds2(nn,1) = bID(nn-1)+1;
                inds2(nn,2) = length(ex);
                dvec(nn) = inds2(nn,2)-inds2(nn,1);
            else
                inds2(nn,1) = bID(nn-1)+1;
                inds2(nn,2) = bID(nn)-1;
                dvec(nn) = inds2(nn,2)-inds2(nn,1);
            end
        end
        [~,II] = sort(dvec,'descend'); % sort the contour groups with longest continuous contours first
        inds2 = inds2(II,:); %Grab the top 5 and ditch the rest
        inds2(6:length(inds2),:) = [];
    else
        inds2 = [1 length(ex)];
    end
    
    ind2use = 1;
    block = inds2(ind2use,1):inds2(ind2use,2);
%     plot(ex(block),ey(block),'k') % To confirm
    
    ex = ex(block);
    ey = ey(block);
    
    
    %% Find where Transect ends intersect with two transects
    % Offshore boundary intersection with transect UTM coordinates
    P2 = InterX([ex;ey],[O(1).x';O(1).y']); % Intersection with first transect
    P1 = InterX([ex;ey],[O(length(O)).x';O(length(O)).y']); % Inersection with last transect
    
    % Find where that intersection is and subsample offshore transect
    % [t1.x,t1.y] = createTransect(O(1).x(1),O(1).y(1),O(1).x(2),O(1).y(2),.2);
    % [t2.x,t2.y] = createTransect(O(length(O)).x(1),O(length(O)).y(1),O(length(O)).x(2),O(length(O)).y(2),.2);
    
    % Find those coordinates in the alongshore transect
    dist = sqrt((P1(1)-ex).^2 + (P1(2)-ey).^2);
    [~, I1] = min(dist);
    dist = sqrt((P2(1)-ex).^2 + (P2(2)-ey).^2);
    [~, I2] = min(dist);
    % Subsample
    ex = ex(I1:I2);
    ey = ey(I1:I2);
    
    %Now make polygon of two cnntours
    px = [T(xx).x' ex T(xx).x(1)];
    py = [T(xx).y' ey T(xx).y(1)];
    
    %% Make flood raster
    floodRaster = SZ - D.Z;
    
    %% Find all points in polygon
    tic
    in = InPolygon(D.X,D.Y,px,py);
    toc
    saveNm = sprintf('%s_inpolygon_SUP_10yr',runName);
    save(saveNm,'in');
    %%
    floodRaster = reshape(floodRaster,size(D.X));
    
end

%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % restoredefaultpath % Necessary for it to work
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % rehash toolboxcache
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % colormap(jet)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % clf
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % IM = load('E:\Abbas\PS_COSMOS\Thesis_Modeling\GeoTiffs\Ruston\ruston_z3.mat');
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % imagesc(IM.xm,IM.ym,IM.im)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % set(gca,'ydir','normal')
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % hold on
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % p = pcolor(D.X,D.Y,floodRaster);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % shading interp
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % caxis([-3 3])
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % c = colorbar;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % axis equal
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % Set to be transparent over land
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % alphaVal = 0.5;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % set(p,'facealpha',alphaVal)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % Get the color data of the object that correponds to the colorbar
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % cdata = c.Face.Texture.CData;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % Change the 4th channel (alpha channel) to 10% of it's initial value (255)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % cdata(end,:) = uint8(alphaVal * cdata(end,:));
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % Ensure that the display respects the alpha channel
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % c.Face.Texture.ColorType = 'truecoloralpha';
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % Update the color data with the new transparency information
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % c.Face.Texture.CData = cdata;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %% Write to KML
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % utmzone = repmat('10 T',length(px),1);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % [lat,lon] = utm2deg(px,py,utmzone);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % Write KML
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % kmlwriteline('swl_SLR280.kml',lat,lon);