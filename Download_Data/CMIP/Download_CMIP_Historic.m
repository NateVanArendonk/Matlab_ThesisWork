clearvars

options = weboptions('Timeout',120); %Set timeout to 60-seconds


% U
fol = 'http://esgdata.gfdl.noaa.gov/thredds/fileServer/gfdl_dataroot/NOAA-GFDL/GFDL-ESM2M/historical/3hr/atmos/3hr/r1i1p1/v20110601/uas/';
for yy = 1860:5:2000
    dstring = sprintf('%d010100-%d123123',yy+1,yy+5);
    fname = sprintf('uas_3hr_GFDL-ESM2M_historical_r1i1p1_%s.nc',dstring);
    
    % Make fpt link
    ftp_file = [fol fname];
    % Save
    save_file = fname;
    
    % Try Once
    mm = 0;
    try
        fprintf('Downloading year %d',yy)
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

% V
fol = 'http://esgdata.gfdl.noaa.gov/thredds/fileServer/gfdl_dataroot/NOAA-GFDL/GFDL-ESM2M/historical/3hr/atmos/3hr/r1i1p1/v20110601/vas/';
for yy = 1860:5:2000
    dstring = sprintf('%d010100-%d123123',yy+1,yy+5);
    fname = sprintf('vas_3hr_GFDL-ESM2M_historical_r1i1p1_%s.nc',dstring);
    
    % Make fpt link
    ftp_file = [fol fname];
    % Save
    save_file = fname;
    
    % Try Once
    mm = 0;
    try
        fprintf('Downloading year %d',yy)
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

% PSFC
fol = 'http://esgdata.gfdl.noaa.gov/thredds/fileServer/gfdl_dataroot/NOAA-GFDL/GFDL-ESM2M/historical/3hr/atmos/3hr/r1i1p1/v20110601/ps/';
for yy = 1860:5:2000
    dstring = sprintf('%d010100-%d123123',yy+1,yy+5);
    fname = sprintf('ps_3hr_GFDL-ESM2M_historical_r1i1p1_%s.nc',dstring);
    
    % Make fpt link
    ftp_file = [fol fname];
    % Save
    save_file = fname;
    
    % Try Once
    mm = 0;
    try
        fprintf('Downloading year %d',yy)
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