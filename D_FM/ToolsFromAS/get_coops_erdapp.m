function data=get_coops_erdapp(stn_id,varargin)
% GET_COOPS_ERDAP - retreive NOAA COOPS data using ERDAPP
%
%   DATA = GET_COOPS_ERDAP(STN_ID) - retrieves water level data for
%       the station specified in STN_ID.  STN_ID should be a character
%       string (eg. '9447130' for Seattle) The user will be prompted for
%       the beginning and ending dates to download and the location to
%       put the resulting file. The data are returned in the structure
%       array, DATA. Times are provided in GMT and water levels in meters.
%
%   OPTIONAL INPUTS
%       Optional inputs are supplied as property pair or as an input
%       structure.
%
%       'service' - type of service requested. Acceptable values are:
%                   'IOOS_SixMin_Verified_Water_Level' (default)
%                   'IOOS_Raw_Water_Level'
%                   'IOOS_Wind'
%       'start_time' - beginning of requested time-series (string)
%       'end_time' - end of requested time-series (string)
%       'datum' - datum of supplied water levels
%       'out_dir' - directory to store the derived .mat file
%
%   EXAMPLES
%       Example 1. User supplies the station id (required) and is prompted
%       for the start and end times, and file destination.
%           stn_id = '9447130'
%           data = get_coops_erdapp(stn_id)
%
%       Example 2. User supplies start and end time through property value
%       pairs.
%           stn_id = '9447130';
%           st = '2015-01-01';
%           et = '2015-01-08';
%           data = get_coops_erdapp(stn_id,'start_time',st,'end_time',et)
%
%       Example 3. User supplies optional inputs as a structure.
%           stn_id = '9447130';
%           opt.start_time = '2015-01-01';
%           opt.end_time = '2015-01-08';
%           opt.datum = 'MLLW';
%           opt.service = 'IOOS_Raw_Water_Level';
%           data = get_coops_erdapp(stn_id,opt);
%
%          %plot the results
%           figure, plot(data.time,data.WL_VALUE)
%           datetick('x')
%
% A map of stations is available from the NOAA Tides and Currents
% page <a href = "https://tidesandcurrents.noaa.gov/">here</a>.
% Information about access to NOAA COOPS data through ERDAPP is
% available <a href = "https://opendap.co-ops.nos.noaa.gov/">here</a>.

% Andrew Stevens, astevens@usgs.gov
% version 1.01
% modified 2/9/2017

p=inputParser;
p.KeepUnmatched=1;

addRequired(p,'stn_id',@ischar);
opts={'service', 'IOOS_SixMin_Verified_Water_Level',{'char'}, {};...
    'start_time', [],      {'char'},    {};...
    'end_time',   [],      {'char'},    {};...
    'datum',     'STND',   {'char'},    {};...
    'out_dir',   [],       {'char'}     {}};

cellfun(@(x)(p.addParameter(x{1},x{2},...
    @(y)(validateattributes(y, x{3},x{4})))),num2cell(opts,2))
p.parse(stn_id,varargin{:});
opt=p.Results;

validatestring(opt.datum,{'MHHW';'MHW';'MTL';'MSL';...
    'MLW';'MLLW';'STND';'NAVD';'IGLD'},'get_coops_erdapp',...
    'datum');
validatestring(opt.service,{'IOOS_SixMin_Verified_Water_Level';...
    'IOOS_Raw_Water_Level';'IOOS_Wind'},'get_coops_erdapp',...
    'service');

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


%only allows downloading 1 month at a time
if dnend-dnstart>30
    dnv=[dnstart:30:dnend,dnend];
    dnstart=dnv(1:end-1);
    dnend=dnv(2:end);
end

%output directory
if isempty(opt.out_dir)
    opt.out_dir = [uigetdir,filesep];
end

if ~exist(opt.out_dir,'dir')
    mkdir(opt.out_dir);
end

%build the URL
url='https://opendap.co-ops.nos.noaa.gov/erddap/tabledap/';
switch opt.service
    case 'IOOS_SixMin_Verified_Water_Level'
        dval='IOOS_SixMin_Verified_Water_Leve';
        url1=sprintf(['?STATION_ID%%2Clatitude%%2Clongitude%%2CDATUM',...
            '%%2Ctime%%2CWL_VALUE&STATION_ID%%3E=%%22',...
            '%s%%22&DATUM%%3E=%%22%s'],opt.stn_id,opt.datum);
    case 'IOOS_Raw_Water_Level'
        dval='IOOS_Raw_Water_Level';
        url1=sprintf(['?STATION_ID%%2Clatitude%%2Clongitude%%2CDATUM',...
            '%%2Ctime%%2CWL_VALUE&STATION_ID%%3E=%%22',...
            '%s%%22&DATUM%%3E=%%22%s'],opt.stn_id,opt.datum);
    case 'IOOS_Wind'
        dval='IOOS_Wind';
        url1=sprintf(['?STATION_ID%%2Clatitude%%2Clongitude',...
            '%%2Ctime%%2CWind_Speed%%2CWind_Direction%%2CWind_Gust',...
            '&STATION_ID%%3E=%%22%s'],opt.stn_id);
end

o = weboptions('CertificateFilename','',...
    'Timeout',30);


gfiles=false(length(dnstart),1);
fname=cell(length(dnstart),1);
for i=1:length(dnstart)
    url2=sprintf(['%%22&BEGIN_DATE%%3E=%%22%s',...
        '%%20%s%%3A%s',...
        '%%22&END_DATE%%3E=%%22%s',...
        '%%20%s%%3A%s%%22'],...
        datestr(dnstart(i),'yyyymmdd'),...
        datestr(dnstart(i),'HH'),datestr(dnstart(i),'MM'),...
        datestr(dnend(i),'yyyymmdd'),...
        datestr(dnend(i),'HH'),datestr(dnend(i),'MM'));
    urlf=[url,opt.service,'.mat',url1,url2];

    fprintf('Downloading file %d of %d, please wait...',i,...
        length(dnstart));
    try
        fname{i}=[opt.stn_id,'_',opt.service,...
            datestr(dnstart(i),'yyyymmdd'),'.mat'];
        websave([opt.out_dir,fname{i}],urlf,o);
        gfiles(i)=true;
        fprintf('Complete!\n'); 
    catch
        fprintf('Problems Downloading File\n'); 
    end
end

if any(gfiles)
    %load the file(s) 
    fnames=fname(gfiles);
    data=arrayfun(@(x)(x.(dval)),...
        cellfun(@(x)(load([opt.out_dir,x])),fnames));
    %concatenate files if necessary
    if length(fnames)>1
        fields=fieldnames(data);
        dc=struct2cell(data)';
        dc2=cellfun(@(x)(cell2mat(dc(:,x))),...
            num2cell(1:length(fields)),'un',0);
        datar=cell2struct(dc2',fields);
        
        %make sure times are unique
        [~,tidx]=unique(datar.time);
        data=structfun(@(x)(x(tidx,:)),datar,'un',0);
    end
    
    %opendap to datenum
    data.time=(data.time./86400) + 719529;
else
    data=[];
  
end



%--subfunction------------------------------------
function out = uigetdate(varargin)
% UIGETDATE  date selection dialog box
%    T = UIGETDATE(D) displays a dialog box in form of a calendar 
%    
%    UIGETDATE expects serial date number or standard MATLAB Date 
%    format (see DATESTR) as input data und returns serial date number 
%    for the selected date and time.
%
%    UIGETDATE by itself uses the current date and time as input data
%
% Example:
%         t = datestr( uigetdate('16-Aug-1974 03:00') )
% 
% See also datevec, datestr, datenum

%   version: v1.0
%   author:  Elmar Tarajan [MCommander@gmx.de]

if nargin == 0
   varargin{1} = now;
end% if

if ~ishandle(varargin{1})
   %
   datvec = datevec(varargin{1});
   %
   scr = get(0,'ScreenSize');
   h.units = 'pixels';
   h.parent = figure(h,'menubar','none', ...
            'numbertitle','off', ...
            'resize','off', ...
            'handlevisibility','on', ...
            'visible','off', ...            
            'WindowStyle','modal', ...
            'Tag','uigetdate', ...
            'position',[ (scr(3:4)- [197 199])/2 197 199 ]);
   %
   pos = [5 5 0 0];
   uicontrol(h,'style','edit','position',pos+[0 0 104 26])
   uicontrol('style','slider','units','pixels','position',pos+[3 2 100 20], ...
             'sliderstep',[.0005 .0005],'min',-10,'max',10,'value',0, ...
             'callback','uigetdate(gcbo,''time'')','UserData',0)
   %
   h.style           = 'edit';
   h.fontweight      = 'bold';
   h.foregroundcolor = [.2 .2 .2];
   uicontrol(h,'enable','inactive','position',pos+[ 17 2 73 20],'Tag','time', ...
               'String',sprintf('%02d:%02d',datvec(4:5)))
   %
   % textbanners
   tmp = [2 20 101 4 ; 2 2 101 3 ; 2 2 3 22 ; 17 2 2 22 ; 88 2 2 22 ; 101 2 2 22 ];
   for i=1:6 ; uicontrol(h,'style','text','position',pos+tmp(i,:)) ; end% for
   %
   uicontrol(h,'style','edit','position',pos+[105 0 84 26],'visible','on')   
   uicontrol(h,'style','pushbutton','position',pos+[108 2 78 21],'Tag','ok', ...
               'CData',repmat(repmat([0.3:0.01:1 1:-0.01:0.3],18,1),[1 1 3]), ...
               'string','ok','Callback','uigetdate(gcbo,''ok'')')
   %
   pos = [5 32 0 0];
   uicontrol(h,'style','edit','position',pos+[0 0 189 136],'enable','inactive','Tag','cday', ...
      'UserData',datvec(3))   
   h.style      = 'pushbutton';
   h.fontweight = 'normal';
   for i=95:-19:0
      for j=0:26:156
         uicontrol(h,'position',pos+[j+3 i+2 27 20],'Enable','off', ...
                     'foregroundcolor',[.2 .2 .2],'Tag','day', ...
                     'callback','uigetdate(gcbo,''day'')')
      end% for
   end% for
   %
   tmp = {'Mon' 'Tue' 'Wed' 'Thu' 'Fri' 'Sat' 'Sun'};
   for j=0:6
      uicontrol(h,'style','text','position',pos+[j*26+4 119 25 13],'string',tmp{j+1}, ...
                  'backgroundcolor',[0.4 0.4 0.4],'foregroundcolor',[.9 .9 .9])         
   end% for
   %
   pos = [5 169 0 0];
   uicontrol(h,'style','edit','position',pos+[0 0 189 26])
   h.style = 'slider';
   uicontrol(h,'position',pos+[3 2 100 20],'sliderstep',[0.00025 1], ...
               'min',-2000,'max',2000,'Value',datvec(2), ...
               'callback','uigetdate(gcbo,''months'')')
   uicontrol(h,'position',pos+[112 2 74 20],'sliderstep',[0.00025 1], ...
               'min',0,'max',4000,'value',datvec(1), ...
               'callback','uigetdate(gcbo,''year'')')
   %
   h.style           = 'edit';
   h.enable          = 'inactive';
   h.fontweight      = 'bold';
   h.foregroundcolor = [.2 .2 .2];
   tmp = {'Januar' 'Februar' 'March' 'April' 'May' 'Juni' 'Juli' ...
          'August' 'September' 'October' 'November' 'December'};
   uicontrol(h,'position',pos+[ 17 2 73 20],'Tag','months','String',tmp{datvec(2)},'Userdata',tmp)
   uicontrol(h,'position',pos+[126 2 47 20],'Tag','year','String',num2str(datvec(1)))
   %
   % textbanners
   h.style = 'text';
   tmp = [2 20 185 4 ; 2 2 185 3 ; 2 2 3 22 ; 17 2 2 22 ; 88 2 2 22 ; ...
      101 2 13 22 ; 126 2 2 22 ; 171 2 2 22 ; 184 2 3 22];
   for i=1:9
      uicontrol(h,'position',pos+tmp(i,:))
   end% for
   %
   set(h.parent,'visible','on')
   setday(varargin{1})
   %
   set(findobj(gcf,'string',num2str(datvec(3))),'CData',geticon)
   %
   uiwait
   try
      out = datenum([num2str( ...
               get(findobj(gcf,'Tag','cday'),'UserData')) '-' ...
               get(findobj(gcf,'Tag','months'),'String') '-' ...
               get(findobj(gcf,'Tag','year'),'String') ' ' ...
               get(findobj(gcf,'Tag','time'),'String') ':00']);
      delete(findobj(0,'Tag','uigetdate'))                       
   catch
      out = [];
      closereq
   end% try
   
   return
end% if

switch varargin{2}
   case 'months'
      h = findobj(gcbf,'Tag','months');
      months = get(h,'UserData');
      set(h,'String',months{mod(get(gcbo,'Value')-1,12)+1})
      set(findobj(gcbf,'Tag','ok'),'Enable','off')      
      %
   case 'year'
      set(findobj(gcbf,'Tag','year'),'String',get(gcbo,'Value'))
      set(findobj(gcbf,'Tag','ok'),'Enable','off')
      %
   case 'day'
      h= findobj(gcf,'Tag','day');
      set(h,'CData',[])

      set(varargin{1},'CData',geticon)
      set(findobj(gcbf,'Tag','cday'),'Userdata',get(varargin{1},'String'))
      set(findobj(gcbf,'Tag','ok'),'Enable','on')
      try 
          uicontrol(h(3)) ;
      end% try %#ok
      return
      %
   case 'time'
      try
         if toc<0.1
            step = get(gcbo,'UserData');
            set(gcbo,'UserData',step+1)
            step = floor(step*sign(get(gcbo,'value'))/2);
         else
            set(gcbo,'UserData',1)
            step = sign(get(gcbo,'value'));
            set(gcbo,'value',0)
         end% if
         %
         handles.time = findobj(gcbf,'Tag','time');
         time = sum(sscanf(get(handles.time,'String'),'%d:%d').*[60;1]);
         time = time+step;
         if time<0
            time = 1439;
         elseif time>1439
            time = 0;
         end% if
         time = sprintf('%02.f:%02.f',floor(time/60),(time/60-floor(time/60))*60);
         set(handles.time,'String',time)
         %
         tic
         return
      catch
         tic
      end% try
      drawnow
      %
   case 'ok'
      uiresume
      return
      %
end% switch
setday(['1-' get(findobj(gcbf,'Tag','months'),'String') '-' ...
             get(findobj(gcbf,'Tag','year'),'String')])
  %
  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setday(datvec)
%-------------------------------------------------------------------------------
datvec = datevec(datvec);
datvec(3) = 1;
%
day = [7 1 2 3 4 5 6];
day = day(weekday(datestr(datvec)));
%
monthend = eomday(datvec(1),datvec(2));
%
ind = [zeros(1,42-monthend-day+1) monthend:-1:1 zeros(1,day-1)];
%
enable = repmat({'on'},42,1);
enable(ind==0) = {'off'};
%
count = strrep(strrep(cellstr(num2str(ind')),' 0',''),' ','');
%
h = findobj(gcf,'Tag','day');
set(h,{'String'},count,{'Enable'},enable,'backgroundcolor',[0.7 0.7 0.7],'CData',[])
set(h(ind~=0),'backgroundcolor',[.925 .922 .9002]);
set(h(ind~=0&repmat([1 1 0 0 0 0 0],1,6)),'backgroundcolor',[1 .8 .8])
  %
  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function icon = geticon
%-------------------------------------------------------------------------------
tmp = [0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ;
       0 0 0 0 1 1 1 1 1 1 1 1 1 1 0 0 1 ; ...
       0 0 0 0 0 1 1 1 1 1 1 1 0 0 0 1 1 ; ...
       1 1 0 0 0 0 1 1 1 1 1 0 0 0 1 1 1 ; ...
       1 1 0 0 0 0 0 1 1 0 0 0 0 1 1 1 1 ; ...
       1 1 1 0 0 0 0 0 0 0 0 0 1 1 1 1 1 ; ...
       1 1 1 1 0 0 0 0 0 0 0 1 1 1 1 1 1 ; ...
       1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 1 1 ; ...
       1 1 1 1 0 0 0 0 0 0 0 0 1 1 1 1 1 ; ...
       1 1 1 0 0 0 0 0 0 0 0 0 0 1 1 1 1 ; ...
       1 1 0 0 0 0 0 1 1 0 0 0 0 0 1 1 1 ; ...
       0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 1 1 ; ...
       0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 1 1 ; ...
       0 0 0 1 1 1 1 1 1 1 1 1 0 0 0 0 1 ; ...
       1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 ];
tmp(tmp==1)=NaN;
tmp(tmp==0)=1;
icon(:,:,1) = tmp;
tmp(tmp==1)=0.25;
icon(:,:,2) = tmp;
tmp(tmp==.25)=0;
icon(:,:,3) = tmp;