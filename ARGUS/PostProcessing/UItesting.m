% Code to load in timestacks and then pick run up points on a graph

clearvars
clc

% Add paths just to be safe
addpath C:\Functions_Matlab\CIRN_Github\UAV-Processing-Toolbox-master
addpath C:\Functions_Matlab\CIRN_Github\PIXel-Toolbox-master
addpath C:\Functions_Matlab

% Clear Cache
loadStack('','clearCache')

% Load Stack
stackFolder = 'E:\ARGUS\Whid_ARGUS\TimeStacks\'; % Name of folder that all stacks are located in
stackName = '1541628000.c2.stack.ras'; % Name of stack file we want to load
[C.p, C.epoch, C.MSC, C.data] = loadStack([stackFolder stackName]); size(C.data); % Load in the stack

% Plot Transects from Image?
plotting = 0;
if plotting
    puvs = plot(C.p.U,C.p.V,'or');
end

%% Load in the variable 'r' used to define the pixels of indices, should have all the meta data we want on each camera/transect
load('E:\ARGUS\Whid_ARGUS\Pix_Timestack_Locations\whidbeyPIX.c2.r.mat')

dinds = find(diff(r.cams.XYZ(:,1))<-5);  % Find each individual transect by finding breaks in vectors

plotting = 0; % Plot to show transects and the extracted transects using dinds above
if plotting
    figure; plot(r.cams.XYZ(:,1),r.cams.XYZ(:,2),'or')
    %figure; plot(diff(r.cams.XYZ(:,1)))
    hold on
    for ii = 1:length(dinds)+1
        if ii == 1
            plot(r.cams.XYZ(1:dinds(ii),1),r.cams.XYZ(1:dinds(ii),2),'.k')
        elseif ii == length(dinds)+1
            plot(r.cams.XYZ(dinds(ii-1)+1:end,1),r.cams.XYZ(dinds(ii-1)+1:end,2),'.k')
        else
            plot(r.cams.XYZ(dinds(ii-1)+1:dinds(ii),1),r.cams.XYZ(dinds(ii-1)+1:dinds(ii),2),'.k')
        end
    end   % Looks ok...
    title('X & Y Points')
end

if plotting
    % Confirm that these indices also hold true for the U and V coordinates
    figure; plot(r.cams.U,r.cams.V,'or')
    %figure; plot(diff(r.cams.XYZ(:,1)))
    hold on
    for ii = 1:length(dinds)+1
        if ii == 1
            plot(r.cams.U(1:dinds(ii)),r.cams.V(1:dinds(ii)),'.k')
        elseif ii == length(dinds)+1
            plot(r.cams.U(dinds(ii-1)+1:end),r.cams.V(dinds(ii-1)+1:end),'.k')
        else
            plot(r.cams.U(dinds(ii-1)+1:dinds(ii)),r.cams.V(dinds(ii-1)+1:dinds(ii)),'.k')
        end
    end
    title('U & V Points')
end

%% Grab and cataloge indices of transects from images
clear inds
for ii = 1:length(dinds)+1
    if ii == 1
        inds(ii,:) = [1 dinds(ii)];
    elseif ii == length(dinds)+1
        inds(ii,:) = [dinds(ii-1)+1 length(r.cams.U)];
    else
        inds(ii,:) = [dinds(ii-1)+1 dinds(ii)];
    end
end
% Populate in indices and XY,UV to the C structure
C.r.inds = inds;
C.r.X = r.cams.XYZ(:,1);
C.r.Y = r.cams.XYZ(:,2);
C.r.U = r.cams.U;
C.r.V = r.cams.V;

% Make images of pixel time stacks and allow user to pick R2 point
% E.g. for line 7, ii = 7
for ii = 10% 1:length(C.r.inds)
    %figure; imshow(C.data(:,inds(ii,1):inds(ii,2)))  % In uv coords.
    
    close all
    size(C.data(:,C.r.inds(ii,1):C.r.inds(ii,2))); % Check the size of data
    t = C.epoch-C.epoch(1); % relative time in seconds
    datetime0 = datetime(datetime(datestr(epoch2Matlab(C.epoch(1))),'TimeZone','UTC'),'TimeZone','America/Los_Angeles'); % Grab first time of image in UTC not Epoch
    ltime0 = datenum(datetime0); datestr(ltime0);  % Local time
    x = C.r.X(C.r.inds(ii,1):C.r.inds(ii,2)); % grabs the cross shore variable for the window of points
    fig = figure; P = pcolor(x,t, C.data(:,C.r.inds(ii,1):C.r.inds(ii,2)));  % In uv coords.
    set(fig,'InvertHardCopy','off')
    set(fig,'Color','k')
    set(gca,'XColor','w','yColor','w','Layer','top')
    
    shading interp; colormap(gray)
    axis ij
    caxis([40 240])
    xlabel('x (m)')
    ylabel(['t (s) after ' datestr(ltime0,'dd-mmm-yyyy HH:MM') 'PDT'])
    title(['x-shore transect ' num2str(ii)],'color','w')
    drawnow
    ylim([500 600])
end
return

%% R2 GUI-like Tool

px = [];
py = [];
sel = 100; % Set the selection equal an empty vector
cx = [];
cy = [];

while(1)
    while sel ~= 32 && sel ~= 8 % once it equals 32 it means that it's time to move on
        [tx,ty,sel] = ginput(1); % get x,y of user input
        cx = vertcat(cx,tx);
        cy = vertcat(cy,ty);
        hold on
        if sel == 1
            p = plot(tx,ty,'ro','MarkerFaceColor','r','MarkerSize',5); % Make marker of plot
        elseif sel == 3 % delete the point
            pm1 = setMarkerColor(p,'r',0.0); % make transparent
            p.MarkerSize = 0.000001; % make super small
        end
    end
    cx(length(cx)) = []; % get rid of the last input that occurs when you hit the space bar or backspace
    cy(length(cy)) = [];
    px = vertcat(px,cx(length(cx))); % user has double clicked so add the x,y to the list
    py = vertcat(py,cy(length(cy)));
    cx = [];
    cy = [];
    if sel == 8
        break
    end
    sel = 1;
end


