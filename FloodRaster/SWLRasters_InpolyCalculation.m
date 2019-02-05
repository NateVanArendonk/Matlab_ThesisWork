%% Load in variables of interest
addpath C:\Functions_Matlab
clear all
close all
clc

M = dir('SWL_NTR_Masks\*.mat');

% Load in All of the Alongshore Transects - to be used later
O = kml2struct('RW_OB_Transects.kml');
for ii = 1:length(O)
    [O(ii).x,O(ii).y] = deg2utm(O(ii).Lat,O(ii).Lon);
end
D = load('RustonSubset_Gridded_Coned.mat');
fprintf('Done Loading Data')
%% Load Xbeach Runs
for xx = 2:length(M)
    %% Make Contour of LowTide line - Arbitrary Depth
    S = load(['SWL_NTR_Masks\' M(xx).name]);
    u_inds = strfind(M(xx).name,'_');
    runName = M(xx).name(u_inds(end)+1:length(M(xx).name));
    dotInds = strfind(runName,'.');
    runName = runName(1:dotInds-1);
    
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
    
    %Now make polygon of two cnntours
    px = [S.x1 fliplr(ex) S.x1(1)];
    py = [S.y1 fliplr(ey) S.y1(1)];
    
    %% Find all points in polygon
    tic
    in = InPolygon(D.X,D.Y,px,py);
    toc
    fprintf('In poly complete %d out of %d\n',xx,length(M))
    saveNm = sprintf('%s_inpolygon_SWL_NTR',runName);
    save(saveNm,'in');
    
end

