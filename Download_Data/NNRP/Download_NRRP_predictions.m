% Download UW WRF NRRP predictions at 12km resolution and 6-hourly
% See readme in http://cses.washington.edu/rocinante/WRF
clearvars

% set ftp location
ftp_loc = 'http://cses.washington.edu/rocinante/WRF/NNRP/wrf_12km/6hrly/';

% time vec
time = datenum(1950,6,13):(1/4):datenum(2010,12,31,18,0,0);
options = weboptions('Timeout',60); %Set timeout to 60-seconds

for tt = 1:length(time)

yr = year(time(tt));
mo = month(time(tt));
dd = day(time(tt));
hr = hour(time(tt));

% % Make directories
% for yy = 1950:2010
%     mkdir(num2str(yy))
% end

% set file name
fname = sprintf('wrfoutp_d02_%04d-%02d-%02d_%02d:00:00',yr,mo,dd,hr);
fsave = sprintf('%d/wrfoutp_d02_%04d%02d%02d%02d.nc',yr,yr,mo,dd,hr);

% Try Once
mm = 0;
try
    data = websave(fsave,sprintf('%s%04d/%s',ftp_loc,yr,fname),options);
catch
    % Try Twice
    try
        data = websave(fsave,sprintf('%s%04d/%s',ftp_loc,yr,fname),options);
    catch
        % Record missing and continue
        mm = mm+1;
        time_missing(mm) = time(tt);
        fprintf('!Could not download %s\n',datestr(time(tt)));
    end
end

end

save('missing_time_stamps.mat','time_missing')
