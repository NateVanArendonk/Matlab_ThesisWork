Location = 'OB';
load('MaxTWL_RustonWay.mat');
load('E:\Abbas\WCRP\Tacoma\Tier3\RecurrenceInterval\OB_RW_TWL_Recurrence.mat');
W = load('E:\Abbas\PS_COSMOS\Thesis_Modeling\Hindcast\Ruston_LUT\OwenRuston_wave_hindcast_SLR0.0FT.mat');
inds = 66:85;
inds2del = [];
for ii = 1:length(E)
    if ~ismember(ii,inds)
        inds2del(end+1) = ii;
    end
end
E(inds2del) = [];
T(inds2del) = [];

events = [];
% Grab all events at or above 10 yr twl 
for ii = 1:length(T)
    inds = T(ii).TWL >= T(ii).TWL_RI10;
    events = horzcat(events,W.time(inds));
end
% Subsample events to get a single event, eliminate events taht were
% multisampled meaning many hours from same event were chosen 
events = unique(events);
events = uniquetol(events,.0001); % This tolerance setting is arbitrary, I played around with it until I could find one that seemed to give me only unique ones 

wnddir = W.wnddir;
spd = W.speed;
time = W.time;
wl = W.twl;
inds2del = [];
for ii = 1:length(events)
    ind = time == events(ii);
    E(ii).spd = spd(ind);
    E(ii).wnddir = wnddir(ind);
    E(ii).wl = wl(ind);
    E(ii).time = time(ind);
    E(ii).date = datestr(time(ind));
    if E(ii).wnddir >= 131 & E(ii).wnddir <= 311 % Heading of Ruston Way 
        inds2del(end+1) = ii;
    end
end
E(inds2del) = [];

% Filter for Speed 
inds2del = [];
for ii = 1:length(E)
    if E(ii).spd < 10
        inds2del(end+1) = ii;
    end
end
E(inds2del) = [];

wl_events = [];
for ii = 1:length(E)
    wl_events(end+1) = E(ii).wl;
end
[~,I] = max(wl_events);
event = E(I);

%% Find where time in TWL timeseries is same as event 
tind = find(W.time == event.time);