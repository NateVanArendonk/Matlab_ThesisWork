clearvars
clc

F = dir('WCRP_Hindcast'); % Get all the folders
F(1:2) = [];
% Grab a single time to just load in the time variable
G = dir(['WCRP_Hindcast/' F(1).name '/*.mat']);
temp = load(['WCRP_Hindcast/' F(1).name '/' G(1).name],'time');
time = temp.time;
years = year(time(1)):1:year(time(end)); % Grab all the years
% Find Each indice equal to each year to speed things up
for ii = 1:length(years)
    yr_inds{ii} = year(time) == years(ii);
end

% Grab already completed Runs
O = dir('WCRP_WaveOut/*.mat');
completedRuns = cell(length(0),1);
for oo = 1:length(O)
    ind = strfind(O(oo).name,'_HsOut');
    O(oo).zoneID = O(oo).name(1:ind-1);
    completedRuns{oo} = O(oo).zoneID;
%     test(oo) = string(O(oo).zoneID);
end




% Loop through and grab the yearly max hsig at each point in each circle 
for ff = 1:length(F)
    G = dir(['WCRP_Hindcast/' F(ff).name '/*.mat']);
    
    % Check to see if the hindcast already has been completed 
    
    if ~ismember(F(ff).name,completedRuns)
        if length(G) > 1
            for gg = 1:length(G)
                W = load(['WCRP_Hindcast/' F(ff).name '/' G(gg).name],'px','py','pz','hs','tp','depth');
                % Now Need to calculate yearly max wave heights
                % Initialize a few variables
                hs_max = zeros(length(W.px),1);
                tp_max = hs_max;
                hs_yr = zeros(length(W.px),length(years));
                tp_yr = hs_yr;
                
                % Go through and grab max and yearly max values
                for ii = 1:length(W.px)
                    [hs_max(ii),MI] = max(W.hs(ii,:)); % Grabs outright max
                    tp_max(ii) = W.tp(ii,MI);
                    for yy = 1:length(years) % Grab yearly max
                        [hs_yr(ii,yy),MI] = max(W.hs(ii,yr_inds{yy}));
                        temp_tp = W.tp(ii,yr_inds{yy});
                        tp_yr(ii,yy) = temp_tp(MI);
                    end
                end
                
                % Grab depths of 8 m or more (depth is positive)
                inds = W.depth <= 8; % anything less than 8m of depth get rid of
                % Take average of yearly max
                hs_yr_avg = nanmean(hs_yr,2);
                tp_yr_avg = nanmean(tp_yr,2);
                % Set anything shallower than 8 meters to zero
                hs_yr_avg(inds) = 0;
                tp_yr_avg(inds) = 0;
                hs_max(inds) = 0;
                tp_max(inds) = 0;
                
                % Store that yearly max in a structure
                H(gg).hs_yr_avg = hs_yr_avg;
                H(gg).tp_yr_avg = tp_yr_avg;
                H(gg).x = W.px;
                H(gg).y = W.py;
                H(gg).z = W.pz;
                H(gg).depth = W.depth;
                clear W
            end
        else
            W = load(['WCRP_Hindcast/' F(ff).name '/' G.name],'px','py','pz','hs','tp','depth');
            % Now Need to calculate yearly max wave heights
            % Initialize a few variables
            hs_max = zeros(length(W.px),1);
            tp_max = hs_max;
            hs_yr = zeros(length(W.px),length(years));
            tp_yr = hs_yr;
            
            % Go through and grab max and yearly max values
            for ii = 1:length(W.px)
                [hs_max(ii),MI] = max(W.hs(ii,:));
                tp_max(ii) = W.tp(ii,MI);
                for yy = 1:length(years)
                    [hs_yr(ii,yy),MI] = max(W.hs(ii,yr_inds{yy}));
                    temp_tp = W.tp(ii,yr_inds{yy});
                    tp_yr(ii,yy) = temp_tp(MI);
                end
            end
            
            % Grab depths of 8 m or more (depth is positive)
            inds = W.depth <= 8;
            % Take average of yearly max
            hs_yr_avg = nanmean(hs_yr,2);
            tp_yr_avg = nanmean(tp_yr,2);
            % Set anything shallower than 8 meters to zero
            hs_yr_avg(inds) = 0;
            tp_yr_avg(inds) = 0;
            hs_max(inds) = 0;
            tp_max(inds) = 0;
            H.hs_yr_avg = hs_yr_avg;
            H.tp_yr_avg = tp_yr_avg;
            H.x = W.px;
            H.y = W.py;
            H.z = W.pz;
            H.depth = W.depth;
            clear W
        end
        saveNm = sprintf('%s_HsOut.mat',F(ff).name);
        save(saveNm,'H')
        movefile(saveNm,'WCRP_WaveOut/')
        fprintf('Completed %d out of %d\n',ff,length(F))
        clear H
    end
end


