% This code will calculate a R2 value at each transect given a TWL and wave
% hindcast
addpath C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\RunUp\RunupFormulas

clearvars


% ------------------- Load in R2 Time series
R = load('VanderMeerR2Out/OB1_R2_rough80.mat');
r2 = R.V2;
hs = R.hs;
twl = R.twl;
ntr = R.ntr;
tp = R.tp;
spd = R.spd;
wnddir = R.wnddir;
TWL = twl+r2;
%% Event Metrics - conditions that brought about max R2 values 
% Find number of events
numEvents = uniquetol(time_max,0.0001); % Searches for unique events within a same day window


% Make a new Structure with Max Events
for ii = 1:length(numEvents)
    E(ii).dtnm = numEvents(ii);
    E(ii).event = datestr(numEvents(ii));
    I = find(time_max == numEvents(ii));
    E(ii).r2 = V2_max(I);
    E(ii).hs = max(hs_max(I));
    E(ii).tp = max(tp_max(I));
    E(ii).wl = max(wl_max(I));
    E(ii).spd = max(spd_max(I));
    E(ii).wnddir = wnddir_max(I(1)); % Didn't want to convert to U and V and then take avg.  
    E(ii).I = I; % Save which transects go with which storm
end

%% Scatter Plot to Explore NTR and Waves with R2 based on TWL
clf
[TWL_sorted,I] = sort(TWL,'descend'); % sort TWL 

sz = 25; % Size of marker 
rng = 1:19000;%length(hs); - basically gets you to ~3m and up TWL

% Get certain number of variables and plot such that lowest values go first 
Sntr = ntr(I(rng));
Shs = hs(I(rng));
Stwl = TWL_sorted(rng);
[~,I] = sort(Stwl,'ascend');

% Scatter plot
scatter(Sntr(I),Shs(I),sz,Stwl(I),'filled')
cb = colorbar;
cb.Label.String = 'TWL [m]';
xlabel('Storm Surge [m]')
ylabel('Wave Height [m]')
grid on 
set(gca,'FontSize',14)



printFig(gcf,'Scatter_NTRvsHsig',[15 11],'png',300)
%% Scatter Plot to Explore NTR and R2 Colored based on TWL
clf
[TWL_sorted,I] = sort(TWL,'descend'); % sort TWL 

sz = 25; % Size of marker 
rng = 1:19000;%length(hs); - basically gets you to ~3m and up TWL

% Get certain number of variables and plot such that lowest values go first 
Sntr = ntr(I(rng));
Sr2 = r2(1,(I(rng)));
Stwl = TWL_sorted(rng);
[~,I] = sort(Stwl,'ascend');

% Scatter plot
scatter(Sntr(I),Sr2(I),sz,Stwl(I),'filled')
cb = colorbar;
cb.Label.String = 'TWL [m]';
xlabel('Storm Surge [m]')
ylabel('Runup [m]')
grid on 
set(gca,'FontSize',14)



printFig(gcf,'Scatter_NTRvsR2',[15 11],'png',300)
%% ID Xbeach Runs - Large Wave Event
[TWL_sorted,I] = sort(TWL,'descend'); % sort TWL to be largest
rng = 1:19000;

% Find max wave conditions 
[max_hs,Ih] = max(hs(I(rng)));

fprintf('Hsig: %.4f, SWL: %.4f, TP: %.2f, TWL: %.2f\n',max_hs, twl(I(rng(Ih))), tp(I(rng(Ih))), TWL(I(rng(Ih))))
fprintf('Spd: %.2f, DIR: %.2f\n',spd(I(rng(Ih))),wnddir(I(rng(Ih))))

%% ID Xbeach Runs - Highest TWL but low water level
clc
[TWL_sorted,I] = sort(TWL,'descend'); % sort TWL to be largest
rng = 1:19000;

fprintf('Hsig: %.4f, SWL: %.4f, TP: %.2f, TWL: %.2f\n',hs(I(rng(1))), twl(I(rng(1))), tp(I(rng(1))), TWL(I(rng(1))))
fprintf('Spd: %.2f, DIR: %.2f\n',spd(I(rng(1))),wnddir(I(rng(1))))
