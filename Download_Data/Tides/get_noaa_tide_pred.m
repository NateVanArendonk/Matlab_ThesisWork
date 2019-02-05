function tide = get_noaa_tide_pred(sta_id,time,interval)
% function tide = get_noaa_tide_pred(sta_id,time,interval)
%
% E.g. Input
% interval = 'hilo'; %6,hourly,hilo
% sta_id = 9447130;
% time = datenum(1991,1,1,0,0,0):(1/24):datenum(2003,1,1,0,0,0);

% Start Stop
date_start = time(1);
date_end = time(end);

% Grab in 5-year increments
yr_start = floor(year(date_start)/5)*5;
yr_end = ceil(year(date_end)/5)*5;

% Loop over 5-year blocks
wl = [];
wl_time = [];
for yr = yr_start:5:yr_end
    T = COOPS_get_tide_predictions(sta_id,[datenum(yr,1,1,0,0,0) datenum(yr+4,12,31,23,0,0)],interval);
    if isempty(T.time)
        error('Interval not available for this station')
    end
    wl = cat(1,wl,T.wl');
    wl_time = cat(1,wl_time,T.time');
end
    
% De-dup time values
[~, I] = unique(wl_time,'first');
wl_time = wl_time(I);
wl = wl(I);

% interp to 1-hour
tide = interp1(wl_time,wl,time,'pchip'); %use pchip, spline increases peak heights

% Check PLOT
% clf
% hold on
% plot(wl_time,wl)
% plot(time,tide)
% datetick('x')

end




function out = COOPS_get_tide_predictions(station,dates,interval)
%COOPS_Wind_Get  Gets COOPS Wind Data from web (http://tidesandcurrents.noaa.gov/index.shtml).
%
%   OUT = COOPS_WIND_GET(station,DATES) gets the wind data for STATION from 
%   (or between and including) DATE(S) where STATION is an NOS id
%   number and DATE(S) is/are Matlab datenum(s) or a year or '*' 
%   (get everything available).
%
%   Input:
%     STATION = num (ex.  station = 8760922)
%     DATE(S) = datenum | [datenum1 datenum2] | year | '*';
%     INTERVAL = 'hourly' (default) | or: 'six minute'
%   Output:
%     OUT = structure:
%           station = BUOY
%           longitude = longitude
%           latitude = latitude
%           time = time
%           temp = temperature (in C)
%     
%
% Dave Thompson (dthompson@usgs.gov)
% edited by Ian Miller (immiller@uw.edu) 25 March 2016
%
%%

% All stations are listed here.
%url = 'http://tidesandcurrents.noaa.gov/station_retrieve.shtml?type=Historic%20Tide%20Data&state=All%20Stations&id1=';
out.name=[];
out.station = station;
station = num2str(station);
if length(dates)==1 % Get one year.
   dates = [datenum(dates,1,1) datenum(dates,12,31)];
end
bdate = datestr(dates(1),'yyyymmdd');
edate = datestr(dates(2),'yyyymmdd');

% Check interval.
if ~any(strcmp(interval,{'6','hourly','hilo'}))
    error('invalid interval')
end

out.interval = interval;

unit = 'metric';
out.units = 'metric';

timezone = 'GMT';   %GMT
out.timezone = 'gmt';

%% read data

url = ['https://tidesandcurrents.noaa.gov/api/datagetter?product=predictions&application=NOS.COOPS.TAC.WL&begin_date=',bdate,...
   '&end_date=',edate,'&datum=MLLW&station=',station,'&time_zone=',timezone,'&units=',unit,'&interval=',interval,'&format=json'];

%http://tidesandcurrents.noaa.gov/api/datagetter?product=predictions&application=NOS.COOPS.TAC.WL&begin_date=19960101&end_date=19970101&datum=MLLW&station=9447130&time_zone=GMT&units=metric&interval=6&format=csv
%http://tidesandcurrents.noaa.gov/api/datagetter?product=predictions&application=NOS.COOPS.TAC.WL&begin_date=19960101&end_date=19961231&datum=MLLW&station=9447130&time_zone=GMT&units=metric&interval=6&format=json


sprintf('Loading %s',url)
   
[tmp,status]=urlread(url);


if status==1

    data=parse_json(tmp);

    if isfield(data,'error')
        out.time=[];
        out.temp=[];
        'No data for that time period.  Moving on.'
    
    else
        
        %assign latitude and longitude of station using metadata returned by NOAA
%         out.latitude=data.metadata.lat;
%         out.longitude=data.metadata.lon;
%         out.name=data.metadata.name;

        % parse data strings
        out.wl = NaN(1,length(data.predictions));
        out.time = NaN(1,length(data.predictions));
        for i=1:length(data.predictions)

            tempdat=data.predictions{i};

            dt=datenum(tempdat.t);
            wl=str2num(tempdat.v);
            
            out.time(i)=dt;
            
            
            %store wind speed to structure
            if ~isempty(wl)
                out.wl(i)=wl;
            else
                out.wl(i)=NaN;
            end            
            clear temp dt wl
        end
    end
elseif status==0
     out.time=[];
     out.wl=[];
     %out.dir=[];
     %out.gus=[];
    'Download Failed!  Moving on.'
end

end

function data = parse_json(string)
% DATA = PARSE_JSON(string)
% This function parses a JSON string and returns a cell array with the
% parsed data. JSON objects are converted to structures and JSON arrays are
% converted to cell arrays.

% F. Glineur, 2009
% (inspired by the JSON parser by Joel Feenstra on MATLAB File Exchange
% (http://www.mathworks.com/matlabcentral/fileexchange/20565) but with 
% faster handling of strings)

pos = 1;
len = length(string);
% String delimiters and escape characters are identified beforehand to improve speed
esc = regexp(string, '["\\]'); index_esc = 1; len_esc = length(esc);

if pos <= len
    switch(next_char)
        case '{'
            data = parse_object;
        case '['
            data = parse_array;
        otherwise
            error_pos('Outer level structure must be an object or an array');
    end
end

    function object = parse_object
        parse_char('{');
        object = [];
        if next_char ~= '}'
            while 1
                str = parse_string;
                if isempty(str)
                    error_pos('Name of value at position %d cannot be empty');
                end
                parse_char(':');
                val = parse_value;
                object.(valid_field(str)) = val;
                if next_char == '}'
                    break;
                end
                parse_char(',');
            end
        end
        parse_char('}');
    end

    function object = parse_array
        parse_char('[');
        object = cell(0, 1);
        if next_char ~= ']'
            while 1
                val = parse_value;
                object{end+1} = val;
                if next_char == ']'
                    break;
                end
                parse_char(',');
            end
        end
        parse_char(']');
    end

    function parse_char(c)
        skip_whitespace;
        if pos > len || string(pos) ~= c
            error_pos(sprintf('Expected %c at position %%d', c));
        else
            pos = pos + 1;
            skip_whitespace;
        end
    end

    function c = next_char
        skip_whitespace;
        if pos > len
            c = [];
        else
            c = string(pos);
        end        
    end
    
    function skip_whitespace
        while pos <= len && isspace(string(pos))
            pos = pos + 1;
        end
    end

     function str = parse_string
        if string(pos) ~= '"'
            error_pos('String starting with " expected at position %d');
        else
            pos = pos + 1;
        end
        str = '';
        while pos <= len
            while index_esc <= len_esc && esc(index_esc) < pos 
                index_esc = index_esc + 1;
            end
            if index_esc > len_esc
                str = [str string(pos:end)];
                pos = len + 1;
                break;
            else
                str = [str string(pos:esc(index_esc)-1)];
                pos = esc(index_esc);
            end
            switch string(pos)
                case '"' 
                    pos = pos + 1;
                    return;
                case '\'
                    if pos+1 > len
                        error_pos('End of file reached right after escape character');
                    end
                    pos = pos + 1;
                    switch string(pos)
                        case {'"' '\' '/'}
                            str(end+1) = string(pos);
                            pos = pos + 1;
                        case {'b' 'f' 'n' 'r' 't'}
                            str(end+1) = sprintf(['\' string(pos)]);
                            pos = pos + 1;
                        case 'u'
                            if pos+4 > len
                                error_pos('End of file reached in escaped unicode character');
                            end
                            str(end+1:end+6) = string(pos-1:pos+4);
                            pos = pos + 5;
                    end
                otherwise % should never happen
                    str(end+1) = string(pos);
                    pos = pos + 1;
            end
        end
        error_pos('End of file while expecting end of string');
    end

    function num = parse_number
        [num, one, err, delta] = sscanf(string(pos:min(len,pos+20)), '%f', 1); % TODO : compare with json(pos:end)
        if ~isempty(err)
            error_pos('Error reading number at position %d');
        end
        pos = pos + delta-1;
    end

    function val = parse_value
        switch(string(pos))
            case '"'
                val = parse_string;
                return;
            case '['
                val = parse_array;
                return;
            case '{'
                val = parse_object;
                return;
            case {'-','0','1','2','3','4','5','6','7','8','9'}
                val = parse_number;
                return;
            case 't'
                if pos+3 <= len && strcmpi(string(pos:pos+3), 'true')
                    val = true;
                    pos = pos + 4;
                    return;
                end
            case 'f'
                if pos+4 <= len && strcmpi(string(pos:pos+4), 'false')
                    val = false;
                    pos = pos + 5;
                    return;
                end
            case 'n'
                if pos+3 <= len && strcmpi(string(pos:pos+3), 'null')
                    val = [];
                    pos = pos + 4;
                    return;
                end
        end
        error_pos('Value expected at position %d');
    end

    function error_pos(msg)
        poss = max(min([pos-15 pos-1 pos pos+20],len),1);
        if poss(3) == poss(2)
            poss(3:4) = poss(2)+[0 -1];         % display nothing after
        end
        msg = [sprintf(msg, pos) ' : ... ' string(poss(1):poss(2)) '<error>' string(poss(3):poss(4)) ' ... '];
        ME = MException('JSONparser:invalidFormat', msg);
        throw(ME);
    end

    function str = valid_field(str)   
    % From MATLAB doc: field names must begin with a letter, which may be
    % followed by any combination of letters, digits, and underscores.
    % Invalid characters will be converted to underscores, and the prefix
    % "alpha_" will be added if first character is not a letter.
        if ~isletter(str(1))
            str = ['alpha_' str];
        end
        str(~isletter(str) & ~('0' <= str & str <= '9')) = '_';   
    end

end

