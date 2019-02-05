% Download UW WRF NRRP predictions at 12km resolution and 6-hourly
% See readme in http://cses.washington.edu/rocinante/WRF
clearvars
addpath C:\Functions_Matlab\time
addpath C:\Functions_Matlab

% set ftp location
ftp_loc = 'http://cses.washington.edu/rocinante/WRF/PNNL_NARR_6km';

% time vec
time = datenum(1981,1,1):1:datenum(2015,12,31);
options = weboptions('Timeout',60); %Set timeout to 60-seconds

% % Make directories
for yy = 1981:2015
    mkdir(num2str(yy))
end

for tt = 1:length(time)

yr = year(time(tt));
mo = month(time(tt));
dd = day(time(tt));

% Make fpt link
%http://cses.washington.edu/rocinante/WRF/PNNL_NARR_6km/1981/data.1981-01-01.nc
ftp_file = sprintf('%s/%d/data.%d-%02d-%02d.nc',ftp_loc,yr,yr,mo,dd);

% Save
save_file = sprintf('%d/data.%d-%02d-%02d.nc',yr,yr,mo,dd);
if exist(save_file,'file') ~= 2
    % Try Once
    mm = 0;
    try
        data = websave(save_file,ftp_file,options);
    catch
        % Try Twice
        try
            data = websave(save_file,ftp_file,options);
        catch
            % Record missing and continue
            mm = mm+1;
            time_missing(mm) = time(tt);
            fprintf('!Could not download %s\n',datestr(time(tt)));
        end
    end
end
end
if exist('time_missing')
    save('missing_time_stamps.mat','time_missing')
end