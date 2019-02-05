% Code to load in timestacks and then pick run up points on a graph

clearvars
clc

% Add paths just to be safe
addpath C:\Functions_Matlab\CIRN_Github\UAV-Processing-Toolbox-master
addpath C:\Functions_Matlab\CIRN_Github\PIXel-Toolbox-master
addpath C:\Functions_Matlab

% Clear Cache
loadStack('','clearCache')

% Load Stack ---- *** CHANGE PATHS
stackFolder = 'E:\ARGUS\Whid_ARGUS\TimeStacks\'; % Name of folder that all stacks are located in
stackName = '1541628000.c2.stack.ras'; % Name of stack file we want to load
rasInd = strfind(stackName,'.ras');
stackNameShort = stackName(1:rasInd-1);
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

% ONLY DO IT FOR ONE TRASECT RIGHT NOW
% Note - Left Click for point
% - Right click to delete current point and change it
% - When happy with point press space bar to log position and move on to
% next R2 crest 


% HERE IS WHERE YOU PICK THE TRANSECT THAT YOU WANT 

for ii = 10% 1:length(C.r.inds)
    %figure; imshow(C.data(:,inds(ii,1):inds(ii,2)))  % In uv coords.
    lowerLimit = 0;
    upperLimit = 100;
    px = [];
    py = [];
    
    % plot the time stacks in 100 second windows - Entire stack is ~600
    % seconds equating to 10 minutes 
    for bb = 1:6
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
        ylim([lowerLimit upperLimit])
        lowerLimit = lowerLimit + 100; % increase window of plotting
        upperLimit = lowerLimit + 100;
        
        %Have user pick R2 points
        sel = 1; % Set the selection equal an empty vector
        cx = [];
        cy = [];
        
        while(1)
            while sel ~= 32 && sel ~= 8 % once it equals 32 it means that it's time to move on
                [tx,ty,sel] = ginput(1); % get x,y,button press of user input
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
            if ~isempty(cx) % logic for if user hits back space without hitting a point first such that cx and cy are empty 
                px = vertcat(px,cx(length(cx))); % user has double clicked so add the x,y to the list
                py = vertcat(py,cy(length(cy)));
            end
            cx = [];
            cy = [];
            if sel == 8
                break % if backspace is hit, gtfo of the while loop 
            end
            sel = 1; % reset sel to start the loop again if we didn't break 
        end
    end
end
%% Save the Camera Data and the pixel points 
C.px = px;
C.py = py;

dirName = sprintf('%d_%d_%d_StackOut',year(ltime0),month(ltime0),day(ltime0));
if ~exist(dirName,'dir')
    mkdir(dirName)
end
saveNm = sprintf('%s_R2_Extracted.mat',stackNameShort);
save(saveNm,'C')
movefile(saveNm,dirName)


%% Now back out the U and V component of each R2 point - Working on this 

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % Should already have this but will set x = the cross shore componenet for
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % indices of that transect denoted by ii, i.e ii = 5 means transect 5
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % x = C.r.X(C.r.inds(ii,1):C.r.inds(ii,2));
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % y = C.r.Y(C.r.inds(ii,1):C.r.inds(ii,2));
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % loop through and find cross shore point closest to each point in the px variable 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % crossX = zeros(length(px),1);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % for pp = 1:length(px)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     [I,~] = findnearest(px(pp),x);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     crossX(pp) = I;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % make temp variable to have entire array of indices for said transect 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % tinds = C.r.inds(ii,1):C.r.inds(ii,2);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % inds = tinds(crossX);





