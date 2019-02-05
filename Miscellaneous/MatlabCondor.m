% This script will try to mimic Condor 
tic
% May need to play around with this section below, not sure if the removal
% of specific cells is univeral 
run_fol = 'E:\Abbas\WCRP\Tacoma\Tier3\XBeach\OB_RW\10yearWL\Contemporary\'; % Folder with different runs housed in individual folders 
run_list = strcat(run_fol,'*');
run_list = dir(run_list);
run_list(1:2) = []; % Get rid of weird shit at beginning 
run_list = nestedSortStruct(run_list,{'date'}); % should sort by date, but if you've been messing around with it it will be all in a differnt order
run_folders = []; 
% Get rid of any non-run files from list of files 
for ii = 1:length(run_list)
    if run_list(ii).isdir == 1
        run_folders(end+1) = ii;
    end
end
run_list = run_list(run_folders);

%% Conodr Equivalent - run multiple XBeach runs and load up next runs once some finish
numQueue = 14;
queue = 1:1:numQueue;
baseDir = pwd;

% Start Each model in the queue
for rr = 1:length(queue)
    cd(strcat(run_list(rr).folder,'/',run_list(rr).name))
    !start xbeach.exe
    cd(baseDir)
end
% Delete queued runs from list 
run_list(queue) = [];

% Act like condor and load up new runs when current ones end 
while ~isempty(run_list)
    % Get number of xbeach windows currently running  and pipe to a text file
    [~,result] = system('tasklist /FI "imagename eq xbeach.exe" /fo table /nh');
    
    % Find how many are running 
    numRunning = length(strfind(result,'.exe')); 
    
    % If runs need to be added, add them as long as we still have plenty of
    % runs to do 
    if numRunning < length(queue) && length(run_list) >= length(queue)
        runs2add = length(queue) - numRunning;
        fprintf('----------------------------------------------------------\n')
        fprintf('She needs %d More RUNS Ma Boy!!!!\n',runs2add)
        fprintf('Time to feed the machine...\n')
        runInds = 1:1:runs2add;
        for ii = 1:length(runInds)
            cd(strcat(run_list(runInds(ii)).folder,'\',run_list(runInds(ii)).name))
            !start xbeach.exe
            cd(baseDir)
            fprintf('Run %s Added\n', run_list(runInds(ii)).name)
        end
        run_list(runInds) = [];
        fprintf('----------------------------------------------------------\n')
    elseif numRunning < length(queue) && length(run_list) < length(queue) % If there are less runs available than the queue length we need to be smart about it 
        slotsAvailable = length(queue) - numRunning;
        runsLeft = length(run_list);
        if runsLeft > slotsAvailable % if there aren't enough spots open, add what we can 
            runs2add = slotsAvailable;
        else % Otherwise add them all 
            runs2add = length(run_list);
        end
        fprintf('----------------------------------------------------------\n')
        fprintf('She needs %d More RUNS Ma Boy!!!!\n',runs2add)
        fprintf('Time to feed the machine...\n')
        runInds = 1:1:runs2add;
        for ii = 1:length(runInds)
            cd(strcat(run_list(runInds(ii)).folder,'\',run_list(runInds(ii)).name))
            !start xbeach.exe
            cd(baseDir)
            fprintf('Run %s Added\n', run_list(runInds(ii)).name)
        end
        run_list(runInds) = [];
        fprintf('----------------------------------------------------------\n')
    end
%     pause(2)
end
fprintf('All Runs Complete\n')
toc