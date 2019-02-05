%% Download Noaa COOPS Tides Verified
clearvars
stn_id = '9443090'; % station number  https://tidesandcurrents.noaa.gov/map/index.shtml?region=Washington
stn_nm = 'Neah';

%'service' - type of service requested. Acceptable values are:
%                   'IOOS_SixMin_Verified_Water_Level' (default)
%                   'IOOS_Raw_Water_Level'
%                   'IOOS_Wind'

%opt.service = 'IOOS_SixMin_Verified_Water_Level';
opt.service =  'IOOS_Barometric_Pressure'; % service chosen from list above 
opt.start_time = '2006-07-11'; % start time
opt.end_time = '2017-09-26'; % end time
opt.datum = 'MLLW'; % Datum, leave alone, will correct below for NAVD88


tides = get_coops_erdapp(stn_id, opt);

%% Convert from MLLW to NAVD88
%%%%%%%%% datum_ example %%%%%%%%%%%%%%%
% [meters] Navd88 = MLLW - mllw2NAVD88
% MLLW: 
%       Seattle MLLW: 2.419 % plug in to vdatum if it does not exist
%       

% EXAMPLE datum_val = 2.419 - 1.7047
%       Seattle = 0.7143 ***Conversion Factor
% Datum conversion page -- https://vdatum.noaa.gov/vdatumweb/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
mllw2navd88 = 0.7150;

tides.WL_VALUE = tides.WL_VALUE - mllw2navd88;  % Tides in MLLW to NAVD88

tides.station_id = stn_id;

tides.DATUM = 'NAVD88';
%% Save file

if strcmp(opt.service,'IOOS_SixMin_Verified_Water_Level')
    file_nm = sprintf('%s_6minV.mat', stn_nm);
elseif strcmp(opt.service,'IOOS_Raw_Water_Level')
    file_nm = sprintf('%s_6minRaw.mat', stn_nm);
elseif strcmp(opt.service,'IOOS_Hourly_Height_Verified_Water_Level')
    file_nm = sprintf('%s_hrV.mat', stn_nm);
else
    file_nm = sprintf('%s_winds.mat', stn_nm);
end

save(file_nm, 'tides')


