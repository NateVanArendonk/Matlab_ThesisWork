function [tideInfo] = findtidalstate2(sampletime, varargin);
% Function to figure out what part of the tide
% cycle a sample was taken. Uses the t_tide
% toolbox available at:
% http://www.eos.ubc.ca/~rich/#T_Tide
%---------------------------------------------
% INPUTS: sampletime: matlab datenum when the
% samples were taken.  Right now this only works on a single
% time (input size is 1 X 1).
%
% The tide station to use can be entered in the
% function command line or done interactively:
%
% OPTIONS: 'station': specifies tidal station to
% use in t_xtide.  If not entered you will be
% prompted to choose a station interactively.
% 'coordinates': you can also search for the
% nearest tidal station in the t_xtide database by
% entering longitude, latitude in the command line
% or enter them interactively....
% 'Graphit': shows a tidal curve for the specified
% time before and after the sample time.
%
%--------------OUTPUT Structure-----------------
%
% tideInfo.sampleht:  tidal height in units 
%       supplied by tidal curve (t_tide)
% tideInfo.sampleslope:  positive or negative 
%       slope calculated from points on curve 
%       immediately surrounding sampletime
% tideInfo.samplestate:  tidal state at 
%       sampletime (ebbing, flooding, turning 
%       to HH, LH, LL, or HL
% tideInfo.springneap:  Spring or Neap period 
%       in cycle
% tideInfo.pcntofHH:  percent of the days range 
%       (LL to HH) that the height represents
%
%--------------Examples------------------------
%
% [tideInfo] = findTidalState2(sampletime);
% or
% [tideInfo = findTidalState2(sampletime,'station',...
%               'Seattle');
% or
% [tideInfo] = findTidalState2(sampletime,'coordinates',...
%               -122.3503 ,47.6218)
% or
% [tideInfo] = findTidalState2(sampletime,'Graphit',2)
%
%--------------FUNCTIONS CALLS------------------
% Internal
%    -t_tide toolbox
%    -findnear.m
%    -peakdetect.m
%   
%-----------------------------------------------
%
% Main code written by Josh Logan USGS 08/2006
% Modified to use t_tide by Andrew Stevens, USGS
% 08/2006


tideOut=[];
tidefreq = 0.54;  %decimal day for 12h and 24min (time between high tides)
lunarday = 1+(48/1440); % 24 hr 48 min 
lunarmonth = 29.53059;

%parse inputs
if nargin < 1
error('Must input the sample time!');
end

%timeVec=[round(sampletime)-15:(1800/86400):round(sampletime)+15];
%Changed back to every 5 minutes, and changed time period to 1 lunar month (Josh Logan)
timeVec=[round(sampletime)-(lunarmonth):(5/(24*60)):round(sampletime)+(lunarmonth)];

for i=1:length(varargin)
    vin = varargin{i};
    if isequal(vin,'station')
        tideStn = varargin{i+1};
        tideTmp=t_xtide(tideStn,timeVec,'format','info');
        tideOut=t_xtide(tideStn,timeVec,'format','raw','units','meters');
    end
    if isequal(vin,'coordinates')
        longitude = varargin{i+1};
        latitude = varargin {i+2};
        tideTmp=t_xtide(longitude,latitude,'format','info');
        tideStn=strtok(tideTmp.station,',');
        tideOut=t_xtide(tideStn,timeVec,'format','raw','units','meters');
    end
    if isequal(vin,'GraphIt')
        int = varargin{i+1};
    end
end




%if not specified, get tidal station name from t_xtide
%first choose station 
if isempty(tideOut)
ButtonName=questdlg('How to find tidal predition station?', ...
                       'Choose prediction station', ...
                       'From List','Latitude/ Longitude','From List');

if strcmp(ButtonName,'From List')==1
        %this is copied straight out of t_xtide to
        %get the tide database
       if ~exist('t_xtide.mat','file'), % Read the harmonics file and make a mat file
  
             fprintf('\n********Can''t find mat-file t_xtide.mat ********\n\n');
             fprintf('Latest version available from http://bel-marduk.unh.edu/xtide/files.html\n\n');   
    
       else
             load t_xtide
       end;
        
       [s,v] = listdlg('PromptString','Select a station:',...
                      'SelectionMode','single',...
                      'ListString',xharm.station);
        tideStn=char(xharm.station(s,:));
        tideStn=strtok(tideStn,',')
        tideOut=t_xtide(tideStn,timeVec,'format','raw','units','meters');
        tideTmp=t_xtide(tideStn,timeVec,'format','info')
else
       prompt  = {'Latitude (N) :','Longitude (E) :'};

            def     = {'0','0'};
            titler   = 'Enter Latitude/Longitude';
            lineNo  = 1;
            dlgresult  = inputdlg(prompt,titler,lineNo,def);
            latitude=str2num(dlgresult{1});
            longitude=str2num(dlgresult{2});
            tideTmp=t_xtide(longitude,latitude,'format','info');
            tideStn=strtok(tideTmp.station,',');
            tideOut=t_xtide(tideStn,timeVec,'format','raw','units','meters');
            
end
end
                   
%Load excel file
mlabdatetime=timeVec;
ht=tideOut';


%Make sure that sampletime is within time period of tidal curve, and that
%tidal curve includes 1/2 lunar day before or after sample time (to
%correctly find HH tide).  i.e. if it's too close to the beginning or end
%of provided tidal curve, then return an error.
if (sampletime < (mlabdatetime(1)+(lunarday/2))) | (sampletime > (mlabdatetime(length(mlabdatetime))-(lunarday/2)))
    sampleht = -9999;
    sampleslope = -9999;
    samplestate = 'ERROR';
    pcntofHH = -9999;
    return
end

%Find indices for all high and low tides using peakdetect
[highidx, lowidx] = peakdetect(ht);
hitimes=mlabdatetime(highidx);
lowtimes=mlabdatetime(lowidx);
hihts = ht(highidx);
lowhts = ht(lowidx);

%Classify all highs as HH (high high), LH (low high) or EH, Equal highs
%Need to create code for cases where LH and HH don't alternate (where there
%are three progressively higher high tides -- 8/28/06 -- quick fix beginning on line 180 )
hitype=cell((length(hihts)-1),1);
for i = 1:(length(hihts)-1)
    if hihts(i) > hihts(i+1)
        hitype(i) = {'HH'};
    elseif hihts(i) < hihts(i+1)
        hitype(i) = {'LH'};
    else 
        hitype(i) = {'EH'}; %Highs are equal
    end
end

%Classify all lows in the same way: LL,HL,EL
lowtype=cell((length(lowhts)-1),1);
for i = 1:(length(lowhts)-1)
    if lowhts(i) > lowhts(i+1)
        lowtype(i) = {'HL'};
    elseif lowhts(i) < lowhts(i+1)
        lowtype(i) = {'LL'};
    else 
        lowtype(i) = {'EH'}; %Highs are equal
    end
end
    
%Find point nearest time on tidal curve
[closest, closestidx] = findnear(mlabdatetime,sampletime,1);

%Get tidal height at sample
sampleht = ht(closestidx);

%Get slope of curve immediately surrounding point
neartimes = mlabdatetime(closestidx-1:closestidx+1);
nearht = ht(closestidx-1:closestidx+1);
dxtime = (neartimes(1)-neartimes(3));
dyht = (nearht(1)-nearht(3));
sampleslope = dyht/dxtime;

%Establish type (ebb, flood, turning)
samplestate = 'undefined';

%if it's exactly on, or within 1 time step of a high or low peak, it's a turning tide
hixx = find((highidx > (closestidx-1)) & (highidx < (closestidx+1)));
if ~isempty(hixx) %means that the closestidx was within 1 time step of a high tide peak
    samplestate = ['turning', char(hitype(hixx))];
end
lowxx = find((lowidx > (closestidx-1)) & (lowidx < (closestidx+1)));
if ~isempty(lowxx) %means that the closestidx was within 1 time step of a low tide peak
    samplestate = ['turning', char(lowtype(lowxx))];
end

%if it isn't exactly on, find direction from slope and find tide type of
%next tide.  Make sure that state hasn't been defined already (above)
if (sampleslope < 0) & (strcmp(samplestate, 'undefined')) == 1  %ebbing, look in low tide peaks
    [peak, peakidx] = findnear(lowidx,closestidx,2);
        if peak < closestidx
            peakidx = peakidx + 1;
        end
     samplestate = ['ebb', char(lowtype(peakidx))];
elseif (sampleslope > 0) & (strcmp(samplestate, 'undefined')) == 1  %flooding, look in high tide peaks
    [peak, peakidx] = findnear(highidx,closestidx,2);
        if peak < closestidx
            peakidx = peakidx + 1;
        end
     samplestate = ['flood', char(hitype(peakidx))];
elseif (strcmp(samplestate, 'undefined')) == 1 %If still undefined
    samplestate = 'error'
end

%Find two high tides surrounding sample time
tmphiidx = find((hitimes > (sampletime-(lunarday/2))) & (hitimes < (sampletime+(lunarday/2))));
if length(tmphiidx) > 2  %Then reduce time by 1 minute on each side of window until you just find two high tides
    minute = datenum([0 0 0 0 1 0]);
    while length(tmphiidx) > 2
        tmphiidx = find((hitimes > (sampletime-(lunarday/2)) + minute) & (hitimes < (sampletime+(lunarday/2)) - minute));
        minute = minute + datenum([0 0 0 0 1 0]);
    end
end
if length(tmphiidx) < 2  %Then increase time by 1 minute on each side of window until you find just two high tides
    minute = datenum([0 0 0 0 1 0]);
    while length(tmphiidx) < 2
        tmphiidx = find((hitimes > (sampletime-(lunarday/2)) - minute) & (hitimes < (sampletime+(lunarday/2)) + minute));
        minute = minute + datenum([0 0 0 0 1 0]);
    end
end

%Find the LL tide for the same time period
tmplowidx = find((lowtimes > (sampletime-(lunarday/2))) & (lowtimes < (sampletime+(lunarday/2))));
if length(tmplowidx) > 2  %Then reduce time by 1 minute on each side of window until you just find two low tides
    minute = datenum([0 0 0 0 1 0]);
    while length(tmplowidx) > 2
        tmplowidx = find((lowtimes > (sampletime-(lunarday/2)) + minute) & (lowtimes < (sampletime+(lunarday/2)) - minute));
        minute = minute + datenum([0 0 0 0 1 0]);
    end
end
if length(tmplowidx) < 2  %Then increase time by 1 minute on each side of window until you find just two low tides
    minute = datenum([0 0 0 0 1 0]);
    while length(tmplowidx) < 2
        tmplowidx = find((lowtimes > (sampletime-(lunarday/2)) - minute) & (lowtimes < (sampletime+(lunarday/2)) + minute));
        minute = minute + datenum([0 0 0 0 1 0]);
    end
end

%Now that two high tides have been found, find the range between HH and LL and calculate the
%percentage of the range that the sampleht represents (above LL).
tmphitypes = hitype(tmphiidx);
tmplowtypes = lowtype(tmplowidx);
HHidx = strmatch('HH',tmphitypes,'exact');
LLidx = strmatch('LL',tmplowtypes,'exact');
tmphihts = hihts(tmphiidx);
tmplowhts = lowhts(tmplowidx);
HHht = tmphihts(HHidx);
if length(HHht)~= 1
    HHht=max(hihts(tmphiidx));
end
LLht = tmplowhts(LLidx);
if length(LLht)~= 1
    LLht=min(lowhts(tmplowidx));
end

%If two HH are found next to each other (in cases where HH and LH don't
%alternate; when high tides get progressively higher for two or more
%cycles) then choose the highest.
if length(HHht) > 1
    HHht = max(HHht);
end
%If two LL are found next to each other (in cases where HL and LL don't
%alternate; when low tides get progressively lower for two or more
%cycles) then choose the lowest.
if length(LLht) > 1
    LLht = min(LLht);
end

daysrange= HHht-LLht;
htfromLL = sampleht-LLht;
pcntofHH = (htfromLL/daysrange) * 100;

%If tidal curve is long enough to calculate a mean range for the month,
%then do it to find whether the period is during spring or neap tides, if
%not, then set springneap variable to 'notcalculated'.
starttime = min(mlabdatetime);
endtime = max(mlabdatetime);

if (endtime - starttime) < 29.53059
    springneap = 'notcalculated';
else
    %get the mean range between HH and LL for the whole time period 
    %*****(this should probably just be for the month around the sample)
    allHHhts = hihts(strmatch('HH',hitype,'exact'));
    allLLhts = lowhts(strmatch('LL',lowtype,'exact'));
    %Get the same lengthe for HH and LL ht vectors
    maxlength = min([length(allHHhts) length(allLLhts)]);
    meanrange = mean(allHHhts(1:maxlength) - allLLhts(1:maxlength));
    
    if daysrange >= meanrange
        springneap = 'spring';
    else
        springneap = 'neap';
    end
end

tideInfo.stationUsed=tideStn;
tideInfo.sampleHeight=sampleht;
tideInfo.tideTimeZone=tideTmp.timezone;
tideInfo.sampleSlope=sampleslope;
tideInfo.sampleState=samplestate;
tideInfo.springNeap=springneap;
tideInfo.pcntofHH=pcntofHH;


if any(strcmp(varargin,'GraphIt'))==1;
    minX=floor(sampletime-(int/2));
    maxX=ceil(sampletime+(int/2));
    fint=round(int)./4;
    
    figure
    plot(timeVec,tideOut)
    hold on
    ylims=get(gca,'ylim');
    %c1=line([sampletime sampletime],[ylims(1) ylims(2)]);
    %set(c1,'linestyle','--','color','k')
    
    %Added by Josh for debugging
    minX = sampletime - lunarmonth/2;
    maxX = sampletime + lunarmonth/2;
%     minX = 732772.006960359 - 1;
%     maxX = 732775.927128183 + 1;
    %----------
    scatter(sampletime,tideInfo.sampleHeight,'*k');
    set(gca,'xlim',[minX maxX],'xtick',[minX:12/24:maxX])
    ylabel('Tide Height (m)','fontsize',12)
    title(['Tide prediction for ',tideStn],'fontsize',12) 
    datetick('x','keeplimits','keepticks')
end
  
