% Code to grab tide predictions using the t_tides toolbox


%staion list http://www.flaterco.com/xtide/locations.html
% datum FAQ - http://www.flaterco.com/xtide/faq.html

%%%%%%%%% datum_val example %%%%%%%%%%%%%%%
% [meters] Navd88 = MLLW - mllw2NAVD88
% MLLW: 
%       Seattle: 2.419

% NAVD88: 
%       Seattle: 1.7047

% EXAMPLE datum_val = 2.419 - 1.7047
%       Seattle = 0.7143 ***Conversion Factor
% Datum conversion page -- https://vdatum.noaa.gov/vdatumweb/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

% Add the t_tide toolbox
addpath C:\Users\ahooshmand\Downloads\t_tide_v1.1
clearvars

% Establish Search Parameters
station = 'seattle';

time_st = datenum(2004,07,21,00,00,00);
time_ed = datenum(2017,07,31,23,54,00);

% Create a 6 minute increasing time vector
TIM = time_st:(6/24/60):time_ed;


tide_t = t_xtide('Seattle (2)', TIM, 'format','raw','unit','meters');

t = [TIM; tide_t];
tide = t';

% Convert to NAVD88
tide(:,2) = tide(:,2) - 0.1122;

tides.time = TIM;
tides.wl = tide(:,2);

file_name = sprintf('%s_tide_predictions.mat', station);

save(file_name, '-struct', 'tides')



%% Create tide preditions using t_xtide

% from Sean

%tims = (now-6/24):(1/48):(now+1);
%tides = t_xtide('La Jolla, Scripps Pier, California',tims,'format','raw','units','feet');