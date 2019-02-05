% Get predictions from NOAA for locations without station data
% Estimate NTR for location (avg nearby stations)
% Estimate TWL for location
clearvars

addpath D:\Functions_Matlab

for sta = 2:4

switch sta
    case 1
        sta_name = 'La Conner';
        sta_id = 9448558;
        mllw2navd88 = 0.461; %[m]
        sta_lat = 48+23.5/60;
        sta_lon = 122+29.8/60;
    case 2
        sta_name = 'Bellingham';
        sta_id = 9449211;
        mllw2navd88 = 0.147; %[m]
        sta_lat = 48+44.7/60;
        sta_lon = 122+29.7/60;
    case 3
        sta_name = 'Nisqually'; %DuPont Wharf
        sta_id = 9446828;
        mllw2navd88 = 1.140; %[m]
        sta_lat = 47+7.1/60;
        sta_lon = 122+39.9/60;
    case 4
        sta_name = 'Olympia';
        sta_id = 9446807;
        mllw2navd88 = 1.211; %[m]
        sta_lat = 47+3.6/60;
        sta_lon = 122+54.2/60;
end

interval = 'hilo';

% Load NTR estimates
ntr_sta = {'Seattle','cherry_point','port_angeles','PortTownsend','Tacoma','La_Push'};
for nn = 1:length(ntr_sta)
    N(nn) = load(sprintf('../ProcessTideData/Output/%s_NAVD88.mat',ntr_sta{nn}));
end

% Use Seattle time series 
time = N(1).time;

% Interp NTR estimates on 1-time series and average
ntr_mean = NaN(5,length(time));
for nn = 1:5 % don't use la push
    temp = interp1(N(nn).time,N(nn).ntr,time);
    ntr_mean(nn,:) = temp;
end
ntr_mean = nanmean(ntr_mean,1);


%% Download tide predictions
% Returns tide in [m, MLLW]
tide_pred = get_noaa_tide_pred(sta_id,time,interval);

% Convert to NAVD88
tide_pred = tide_pred - mllw2navd88;

% Estimate TWL
twl = tide_pred + ntr_mean;


%% Save to file

save(sprintf('Synthetic_TWL_%d_NAVD88.mat',sta_id),'time','twl','ntr_mean','tide_pred','sta_name','sta_id','sta_lat','sta_lon');






end





return


%% Quick plot of NTR's to compare
% Note, La Push is on open coast, don't use in average! 
clf

subplot(211)
hold on
for nn=1:6
    plot(N(nn).time,N(nn).ntr,'.-')
    %plot(N(nn).time,N(nn).ntr_simple,'.-')
end
plot(time,ntr_mean,'k')
xlim([datenum(2006,12,10) datenum(2007,1,7)])
datetick('x','mm/dd/yyyy','keeplimits')
legend('Seattle','cherry point','port angeles','PortTownsend','Tacoma','La Push','Location','Northwest')
ylabel('NTR [m]')
set(gcf,'Renderer','painters')

nn = 6;
subplot(212)
hold on
plot(N(nn).time,N(nn).tide,'.-')
plot(N(nn).time,N(nn).tide_pred,'.-')
xlim([datenum(2006,12,10) datenum(2007,1,7)])
datetick('x','mm/dd/yyyy','keeplimits')
ylabel('La Push WL [m, NAVD88]')
legend('Obs','Pred')

printFig(gcf,'NTR_ShortTimeSeries_Comparison',[11 8],'pdf')
