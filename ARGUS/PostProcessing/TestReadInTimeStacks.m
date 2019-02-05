% Test code to read in a timestack image from ARGUS 
clearvars 
clc

% 2. Clear Cache 
loadStack('','clearCache')
return 
% 3. Load Stack 
stackFolder = 'E:\ARGUS\Whid_ARGUS\TimeStacks\';
stackName = '1541628000.c2.stack.ras';
%[P,E,M,D] = loadStack(stackName,'UV'); % Doesn't Work
[P,E,M,D] = loadStack([stackFolder stackName]); % Works 
% P = loadStack(stackName,'info'); % Works 

% Put in Argus time not Epoch time 
P.UTC = datenum(datetime(P.syncTime,'convertfrom','posixtime'));
P.PST = P.UTC - 8/24;

%% Plot
clf
restoredefaultpath % Necessary for it to work
rehash toolboxcache
imagesc(P.U,P.V,D)
title(datestr(P.PST))
colormap(bone)
set(gca,'ydir','normal')
hold on
