function data =  get_usgs_waterdata(site,varargin)
% GET_USGS_WATERDATA - retreive USGS waterdata using web services
%
%   DATA = GET_USGS_WATERDATA(SITE) - retrieves discharge and gauge height
%       data for the station specified in SITE.  SITE should be a character
%       string (eg. '12045500' for Elwha River at McDonald Bridge) 
%       specifying a single station. The user will be prompted for the 
%       beginning and ending dates to download. The data are returned 
%       in the structure array, DATA. Times are provided in GMT, discharge
%       in m3/s and gauge height in meters.
%
%   OPTIONAL INPUTS
%       Optional inputs are supplied as property pair or as an input
%       structure.
%
%       'start_time' - beginning of requested time-series (string)
%       'end_time' - end of requested time-series (string)
%
%   EXAMPLES
%       Example 1. User supplies the station id (required) and is prompted
%       for the start and end times.
%           site = '12045500';
%           data = get_usgs_waterdata(site)
%
%       Example 2. User supplies start and end time through property value
%       pairs. Get last week of data.
%           site = '12045500';
%           st = datestr(now-7);
%           et = datestr(now);
%           data = get_usgs_waterdata(site,'start_time',st,'end_time',et)
%
%       Example 3. User supplies optional inputs as a structure.
%           site = '12045500';
%           opt.start_time = datestr(now-7);
%           opt.end_time = datestr(now);
%           data = get_usgs_waterdata(site,opt);
%
%          %plot the results
%           figure,
%           ax(1)=subplot(211);
%           plot(data.mtime,data.cms)
%           ylabel('Discharge (m^3/s)')
%
%           ax(2)=subplot(212);
%           plot(data.mtime,data.ht_m)
%           ylabel('Gauge Height (m)')
%
%           linkaxes(ax,'x')
%           set(ax(1),'xticklabel',[])
%           datetick('x')
%
%       Example 4. Download the data for the last week for a group of sites
%       in the Puget Sound area
%           sites={'12045500';... %Elwha
%                  '12056500';... %Skokomish
%                  '12089500';... %Nisqually
%                  '12096505';... %Puyallup
%                  '12167000';... %Stilliguamish
%                  '12200500';... %Skagit
%                  '12213100'};   %Nooksack
%           opt.start_time = datestr(now-7);
%           opt.end_time = datestr(now);           
%           data=cellfun(@(x)(get_usgs_waterdata(x,opt)),sites);
%
% More information on USGS water data can be found <a href = "https://waterdata.usgs.gov/usa/nwis/">here</a>.
% Information about access to other USGS water data products through web
% services is available <a href = "https://waterservices.usgs.gov/">here</a>.

% Andrew Stevens, astevens@usgs.gov
% version 1.00
% modified 12/1/2017


p=inputParser;
p.KeepUnmatched=1;

addRequired(p,'site',@ischar);
opts={'start_time', [],      {'char'},               {};...
      'end_time',   [],      {'char'},               {}};

cellfun(@(x)(p.addParameter(x{1},x{2},...
    @(y)(validateattributes(y, x{3},x{4})))),num2cell(opts,2))
p.parse(site,varargin{:});
opt=p.Results;



%start and end times
if isempty(opt.start_time)
    fprintf('Pick a start time.\n');
    dnstart=uigetdate(now);
else
    dnstart=datenum(opt.start_time);
end
if isempty(opt.end_time)
    fprintf('Pick an end time.\n');
    dnend=uigetdate(dnstart);
else
    dnend=datenum(opt.end_time);
end

if dnstart>now
    error('Pick a time period starting in the past.')
end

if dnend-dnstart>365
    dnv=[dnstart:365:dnend,dnend];
    dnstart=dnv(1:end-1);
    dnend=dnv(2:end);
end

%only allows downloading 1 year at a time
mtime=cell(length(dnstart));
cms=cell(length(dnstart));
ht=cell(length(dnstart));


webopts=weboptions('Timeout',20);
for i=1:length(dnstart)
    
    url=['https://nwis.waterservices.usgs.gov/nwis/iv/?format=json',...
        '&sites=',opt.site,...
        '&startDT=',datestr(dnstart(i),'yyyy-mm-dd'),...
        '&endDT=',datestr(dnend(i),'yyyy-mm-dd'),...
        '&parameterCd=00060'...
        '&siteStatus=all'];
    
    rdata=webread(url,webopts);
    
    %only allow one station per request
    if i==1
        if ~isempty(rdata.value.timeSeries)
            sourceInfo=rdata.value.timeSeries(1).sourceInfo;
            data.siteName=sourceInfo.siteName;
            data.siteId=sourceInfo.siteCode.value;
            data.site_lon=sourceInfo.geoLocation.geogLocation.longitude;
            data.site_lt=sourceInfo.geoLocation.geogLocation.latitude;
        else
            data=[];
            return
        end
    end

    cms{i}=arrayfun(@(x)(str2double(x.value).*0.0283168),...
        rdata.value.timeSeries(1).values.value);
    %cfs to metric
    
%     ht{i}=arrayfun(@(x)(str2double(x.value).*0.3048),...
%         rdata.value.timeSeries(2).values.value);
    %ft to m

    
    %convert to datetimes
    ds=arrayfun(@(x)(x.dateTime(1:23)),...
        rdata.value.timeSeries(1).values.value,'un',0);
    tz=arrayfun(@(x)(str2double(x.dateTime(25:26))),...
        rdata.value.timeSeries(1).values.value);
    mtime{i}=datenum(ds,'yyyy-mm-ddTHH:MM:SS.FFF')-(tz/24); %gmt
end

mt=cell2mat(mtime);
dcms=cell2mat(cms);
% data.ht_m=cell2mat(ht);
[data.mtime,midx]=unique(mt);
data.cms=dcms(midx,1);
% data.ht_m=dht(midx,1);




