% Dowload NARR 
% Extract u,v,prtmsl
% Save in .mat file
% Overwrite grib files (6TB of total data downloaded)
% Loc
%   https://nomads.ncdc.noaa.gov/data/narr/
% Lat/Lon
%   https://nomads.ncdc.noaa.gov/data/narr/latlon.g221.txt

clearvars 

addpath E:\NARR\nctoolbox-master
setup_nctoolbox


% set ftp location
ftp_loc = 'https://nomads.ncdc.noaa.gov/data/narr';
options = weboptions('Timeout',60); %Set timeout to 60-seconds

% time vec
time = datenum(1979,1,1):(1/8):datenum(2015,12,31);

for tt = 1121:length(time)

    fprintf('Downloading %s\n',datestr(time(tt)));
    
% Create http url
fol = sprintf('%s/%s',datestr(time(tt),'yyyymm'),datestr(time(tt),'yyyymmdd'));
dstring = datestr(time(tt),'yyyymmdd_HHMM');
http_file = sprintf('%s/%s/narr-a_221_%s_000.grb',ftp_loc,fol,dstring);
save_file = 'C:\Users\ahooshmand\Documents\MATLAB\test.grb';

% Download, keep trying, never stops
fail = true;
count = 0;
while fail
    try 
        count = count+1;
        data = websave(save_file,http_file,options);
        fail = false;
    catch
        fprintf('Failed, attemp %d\n',count)
    end
end

% Load grib
nco = ncgeodataset(save_file);

% Height above ground, 30m, 10m, 2m
myvar = 'height_above_ground';
temp = nco.geovariable(myvar);
hag = temp.data(:); 

% Time
% myvar = 'time';         % 'Hour since 1979-01-01T00:00:00Z'
% temp = nco.geovariable(myvar);
% grib_time = temp.data(:); % Height above ground, 30m, 10m, 2m               
% grib_time_unit = temp.attributes{1,2};
time_step = time(tt);

% Velocity at 10m
myvar = 'u_wind_height_above_ground'; %vel [m/s], select 2, 10-m above ground
temp = nco.geovariable(myvar);
u = squeeze(temp.data(1,2,:,:));

myvar = 'v_wind_height_above_ground'; %vel [m/s], select 2, 10-m above ground
temp = nco.geovariable(myvar);
v = squeeze(temp.data(1,2,:,:));

% PRTMSL
myvar = 'Pressure_reduced_to_MSL_msl'; %Pressure [Pa]
temp = nco.geovariable(myvar);
prtmsl = squeeze(temp.data(1,:,:));

% Make folder
if ~exist(fol,'dir')
    mkdir(fol)
end
    
% Save
save_mat = sprintf('%s/narr-a_221_%s',fol,dstring);
save(save_mat,'u','v','prtmsl','hag','time_step')

delete('C:\Users\ahooshmand\Documents\MATLAB\test.grb.gbx9')
% delete('C:\Users\ahooshmand\Documents\MATLAB\test.grb.ncx')

fprintf('Completed %d out of %d - Moving On...\n',tt,length(time));
end

return
