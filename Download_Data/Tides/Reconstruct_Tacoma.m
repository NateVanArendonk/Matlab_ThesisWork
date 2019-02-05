% Reconstruct water levels from pred and NTR

clearvars

T = load('../NTR/Tacoma_NTR_NAVD88.mat');
S = load('../NTR/seattle_NTR_NAVD88.mat');

% Get long-term Tacoma predictions
% default MLLW, 6-min for predictions
return
tic % like 45 minutes, maybe 1-hour
[P.time,P.pred,P.meta] = getNOAAtide_hourly(datestr(S.time(1),'yyyymmdd'),datestr(S.time(end),'yyyymmdd'),num2str(T.id));
toc
P.pred = P.pred - T.mllw2navd88;

% Interp on same time axis
O.time = S.time;
O.pred = interp1(P.time,P.pred,O.time);
O.twl = O.pred + S.ntr; 

%% Compare reconstructed to obs
clf
twl_recon = interp1(O.time,O.twl,T.time);

clf
hold on
plot(T.time,T.pred)
plot(T.time,T.obs)
plot(T.time,twl_recon)

xlim([datenum(2017,1,1) datenum(2017,2,1)])
legend('pred','obs','reconstructed')

%% Patch obs and recon
I = findnearest(T.time(1),O.time);
O.twl_patched = O.twl;
O.twl_patched(I:end) = interp1(T.time,T.obs,O.time(I:end));

%% Save
O.id = T.id;
O.name = T.name;
O.ntr = S.ntr;

save(sprintf('%s_Reconstruct_NAVD88.mat',O.name),'-struct','O');

%%
clf
hold on
plot(O.time,O.twl)
plot(T.time,T.obs)
legend('Reconstructed','Observed')
grid on
box on
%plot(O.pred)
ylim([3 4])
datetick
ylabel('Water Level [m, NAVD88]')
title('Tacoma')
xlim([datenum(1900,1,1) datenum(2017,1,1)])
printFig(gcf,'TimeSeries_Tacoma_highTWL',[10 6],'png')

