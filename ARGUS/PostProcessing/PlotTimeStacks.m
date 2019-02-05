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
[C.p, C.epoch, C.MSC, C.data] = loadStack([stackFolder stackName]); size(C.data) % Load in the stack

% Plot Transects from Image?
plotting = 0;
if plotting
    puvs = plot(C.p.U,C.p.V,'or');
end

%% Load in the variable 'r' used to define the pixels of indices, should have all the meta data we want on each camera/transect 
load('E:\ARGUS\Whid_ARGUS\Pix_Timestack_Locations\whidbeyPIX.c2.r.mat')

dinds = find(diff(r.cams.XYZ(:,1))<-5);  % Find each individual transect by finding breaks in vectors

plotting = 1; % Plot to show transects and the extracted transects using dinds above
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
end
title('X & Y Points')

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
end
title('U & V Points')
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

% Make images of pixel time stacks 
% E.g. for line 7, ii = 7
for ii = 5% 1:length(C.r.inds)
    %figure; imshow(C.data(:,inds(ii,1):inds(ii,2)))  % In uv coords.
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
    saveNm = ['whidbey.tran' num2str(ii,'%02.0f') '.' num2str(C.epoch(1),'%10.0f') '.C.png'];
    printFig(gcf,saveNm,[4 11],'png',300)
    pause(.1)
    ylim([0 180])
    saveNm = ['whidbey.tran' num2str(ii,'%02.0f') '.' num2str(C.epoch(1),'%10.0f') '.C.ZOOM.png'];
    printFig(gcf,saveNm,[4 11],'png',300)
    close gcf
end


