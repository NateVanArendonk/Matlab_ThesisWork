function pred = download_tidePredict(station_id, dates, datum_val, interval)
% Script to download tide predictions for NOAA COOPS stations 

% station_id = ID of station, ex. Seattle, Wa = 9447130
% dates = dates to download ex [20160101, 20170101]
% ******Note 9 year max for download
% datum_val = Conversion of water level datum to NAVD88
% interval = hourly or six minute: input either -> 'h' or '6'


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
%% Build NOAA URL
% First create a weboptions object and set the Timeout property to 30
options = weboptions;
options.Timeout = Inf;



url1 = ['https://tidesandcurrents.noaa.gov/cgi-bin/predictiondownload.cgi?&stnid=',station_id];
url2 = ['&threshold=&thresholdDirection=&bdate=',num2str(dates(1)),'&edate=',num2str(dates(2))];
url3 = ['&units=metric&timezone=GMT&datum=MLLW&interval=',interval,'&clock=12hour&type=txt&annual=false'];
    
url = strcat(url1, url2, url3);
%file_nm = sprintf('%s_predicted.txt', station_id);

%download_data
web_pg = webread(url, options); % grab web data from URL using weboptions supplied above

if interval == 'h'
    web_data = web_pg(406:end); % get rid of header
else
    web_data = web_pg(411:end); % get ride of header for six minute data

% Parse data
data_parse = get_tokens(web_data, '\t');
% convert data into columns of time, time, day, tide
data = reshape(data_parse, [4, length(data_parse)/4])';

% Grab the size of the data matrix
[row, col] = size(data);

%% Populate variables with data

% First Initialize a few empty variables
count = 0;
tdate = NaN(1,row);
tide = NaN(1,row);


format = 'yyyy/mm/dd HH:MM'; % format for datenum conversion

for r = 1:row
    temp_dt = data(r); % populate a temporary variable with the data, first variable in the data set
    temp_dt = temp_dt{1}; % Convert to string
    for c = 3:col
        if c == 3 % First deal with the date
            tt = data(r,c); % grab the data for hour and minute
            tt = tt{1}; % Convert to string
            date_str = sprintf('%s %s', temp_dt, tt(1:5)); % concatenate and create a time string for that time 
            tdate(r) = datenum(date_str, format); % add datenum to variable
        else % Now deal with the tide value
            bb = data(r,c); % grab tide data
            tide_str = bb{1}; % grab tide data that is a string
            tide(r) = str2num(tide_str) - datum_val; % Add the number to the tide variable and convert to NAVD88
        end
    end   
end
% Create structure of data
    pred.tides = tide';
    pred.time = tdate';
end

