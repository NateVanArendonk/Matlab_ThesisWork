% This script will try to mimic Condor
tic

% Get list of all the run folders 
RunList = 'C:\Users\egrossman\Desktop\XBeach_LUT_sites\XBRuns\';
Runs = dir([RunList '*']);
Runs(1:2) = [];

for zz = 1:length(Runs) % Go through each run folder and run condor on all the runs 
run_fol = [RunList Runs(zz).name '\']; % Folder with different runs housed in individual folders
xbeach_path = 'C:\Users\egrossman\Desktop\XBeach_LUT_sites\';
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
numQueue = 26;
queue = 1:1:numQueue;
baseDir = pwd;

% Start Each model in the queue
for rr = 1:length(queue)
    cd(strcat(run_list(rr).folder,'/',run_list(rr).name))
    copyfile([xbeach_path 'xbeach.exe'],strcat(run_list(rr).folder,'/',run_list(rr).name));
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
            copyfile([xbeach_path 'xbeach.exe'],strcat(run_list(runInds(ii)).folder,'\',run_list(runInds(ii)).name));
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
            copyfile([xbeach_path 'xbeach.exe'],strcat(run_list(runInds(ii)).folder,'\',run_list(runInds(ii)).name));
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
end
return
%% Go in and Remove all of the xbeach.exe files
run_fol = 'C:\Users\egrossman\Desktop\XBeach_LUT_sites\XBRuns\ID10_xyz\'; % Folder with different runs housed in individual folders
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

for ii = 1:length(run_list)
    cd(strcat(run_list(rr).folder,'/',run_list(rr).name))
    delete('xbeach.exe')
    cd(baseDir)
end
