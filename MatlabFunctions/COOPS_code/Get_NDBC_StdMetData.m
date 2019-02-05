close all
clear all

%collate NDBC historical buoy data (only accesses data for the nearest whole year). 
%DOES NOT use NetCDF instead downloads .txt summary files served by NDBC
%see http://www.ndbc.noaa.gov/measdes.shtml for measurement descriptions and units
%script also calculates wave power (wave energy flux) and adds it to the
%output
% Ian Miller, immiller@uw.edu
%VERSION UPDATE:  16 December 2016

stname='46041'; %46088=New Dungeness Buoy; 46041 Cape Elizabeth buoy; desw1 Destruction Island; 46087 Neah Bay
yrstart=1970; %defaults to downloading ALL data but these can be manually changed
temp=datevec(today);
yrend=temp(1); clear temp

%initialize output structure
StdMetData=struct('station',stname,'time',[],'wvht',[],'dpd',[],'apd',[],...
    'mwvd',[],'wavepower',[],'winddir',[],'windspd',[],'gust',[],'pres',[],'airtemp',[],...
    'watertemp',[]);

for i=yrstart:1:yrend
   url=['http://www.ndbc.noaa.gov/view_text_file.php?filename=' stname 'h' num2str(i) '.txt.gz&dir=data/historical/stdmet/'];
   
   sprintf('Reading %s data',num2str(i))
   [inr,status]=urlread(url);
   status
   
   if status==1
   
       %read first line to get variable headers
       checkline=textscan(inr,'%s',1,'Delimiter','/n');

       %get number of variables in file and build appropriate format string for
       %textscan
       vars=strsplit(char(checkline{1}),' ');

       if length(vars)==14
           formstr='%f%f%f%f%f%f%f%f%f%f%f%f%f%f';
       elseif length(vars)==15
           formstr='%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f';
       elseif length(vars)==16
           formstr='%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f';
       elseif length(vars)==17
           formstr='%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f';
       elseif length(vars)==18
           formstr='%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f';
       end

       %read in rest of data file as a series of strings
       data=textscan(inr,char(formstr),'headerlines',2,'EndofLine','\r\n');
       clear inr url %no longer needed so delete

       %% match date/time data in vars to output variable names for date/time

       if ~isempty(find(strncmp('YY',vars,2))) %this takes into account varying YR variable headers in NDBC text files
            yrs=data{find(strncmp('YY',vars,2))}; 
       else
           yrs=data{find(strncmp('#YY',vars,3))};
       end

       if nanmean(yrs)<100
           yrs=yrs+1900;
       elseif nanmean(yrs)<20
           yrs=yrs+2000;
       end

       mos=data{strmatch('MM',vars)}; 
       dy=data{strmatch('DD',vars)}; hrs=data{strmatch('hh',vars)};

       if strmatch('mm',vars) %this accounts for variation in ndbc text files - some include mintues, some don't
           dts=datenum(yrs,mos,dy,hrs,data{strmatch('mm',vars)},0);
       else
           dts=datenum(yrs,mos,dy,hrs,0,0);
       end

       StdMetData.time=[StdMetData.time; dts];

       %% now run through variable and assign to structure where possible
       
       %wave height
       if ~isempty(find(strncmp('WVHT',vars,4)))
           temp=data{find(strncmp('WVHT',vars,4))}; %wave height data
           temp(find(temp==99))=NaN;
           StdMetData.wvht=[StdMetData.wvht; temp];
           clear temp
       else
           StdMetData.wvht=[StdMetData.wvht; NaN(length(dts),1)];;
       end

       %dominant period
       if ~isempty(find(strncmp('DPD',vars,3)))
           temp=data{find(strncmp('DPD',vars,3))}; %wave height data
           temp(find(temp==99))=NaN;
           StdMetData.dpd=[StdMetData.dpd; temp];
           clear temp
       else
           StdMetData.dpd=[StdMetData.dpd; NaN(length(dts),1)];;
       end

       %average period
       if ~isempty(find(strncmp('APD',vars,3)))
           temp=data{find(strncmp('APD',vars,3))}; %wave height data
           temp(find(temp==99))=NaN;
           StdMetData.apd=[StdMetData.apd; temp];
           clear temp
       else
           StdMetData.apd=[StdMetData.apd; NaN(length(dts),1)];;
       end

       %mean wave direction
       if ~isempty(find(strncmp('MWD',vars,3)))
           temp=data{find(strncmp('MWD',vars,3))}; %wave height data
           temp(find(temp==999))=NaN;
           StdMetData.mwvd=[StdMetData.mwvd; temp];
           clear temp
       else
           StdMetData.mwvd=[StdMetData.mwvd; NaN(length(dts),1)];;
       end

       %wind direction
       if ~isempty(find(strncmp('WDIR',vars,4)))
           temp=data{find(strncmp('WDIR',vars,4))}; %wave height data
           temp(find(temp==999))=NaN;
           StdMetData.winddir=[StdMetData.winddir; temp];
           clear temp
       elseif ~isempty(find(strncmp('WD',vars,2)))
           temp=data{find(strncmp('WD',vars',2))};
           temp(find(temp==999))=NaN;
           StdMetData.winddir=[StdMetData.winddir; temp];
           clear temp
       else
           StdMetData.winddir=[StdMetData.winddir; NaN(length(dts),1)];
       end
       
       %wind speed
       if ~isempty(find(strncmp('WSPD',vars,4)))
           temp=data{find(strncmp('WSPD',vars,4))}; %wave height data
           temp(find(temp==99))=NaN;
           StdMetData.windspd=[StdMetData.windspd; temp];
           clear temp
       else
           StdMetData.windspd=[StdMetData.windspd; NaN(length(dts),1)];;
       end
       
       %wind gust
       if ~isempty(find(strncmp('GST',vars,3)))
           temp=data{find(strncmp('GST',vars,3))}; %wave height data
           temp(find(temp==99))=NaN;
           StdMetData.gust=[StdMetData.gust; temp];
           clear temp
       else
           StdMetData.gust=[StdMetData.gust; NaN(length(dts),1)];;
       end
       
       %atmospheric pressure
       if ~isempty(find(strncmp('PRES',vars,4)))
           temp=data{find(strncmp('PRES',vars,4))}; %wave height data
           temp(find(temp==9999))=NaN;
           StdMetData.pres=[StdMetData.pres; temp];
           clear temp
       elseif ~isempty(find(strncmp('BAR',vars,3)))
           temp=data{find(strncmp('BAR',vars',3))};
           temp(find(temp==9999))=NaN;
           StdMetData.pres=[StdMetData.pres; temp];
           clear temp
       else
           StdMetData.pres=[StdMetData.pres; NaN(length(dts),1)];
       end
       
       %air temperature
       if ~isempty(find(strncmp('ATMP',vars,4)))
           temp=data{find(strncmp('ATMP',vars,4))}; %wave height data
           temp(find(temp==999))=NaN;
           StdMetData.airtemp=[StdMetData.airtemp; temp];
           clear temp
       else
           StdMetData.airtemp=[StdMetData.airtemp; NaN(length(dts),1)];;
       end
       
       % water temperature
       if ~isempty(find(strncmp('WTMP',vars,4)))
           temp=data{find(strncmp('WTMP',vars,4))}; %wave height data
           temp(find(temp==999))=NaN;
           StdMetData.watertemp=[StdMetData.watertemp; temp];
           clear temp
       else
           StdMetData.watertemp=[StdMetData.watertemp; NaN(length(dts),1)];;
       end
   elseif status==0
   end
end

%% add more recent data

%monthly data - NEED TO COMPLETE

% for i=1:12
%    url=['http://www.ndbc.noaa.gov/view_text_file.php?filename=' stname 'h' num2str(i) '.txt.gz&dir=data/historical/stdmet/'];
%    http://www.ndbc.noaa.gov/view_text_file.php?filename=4604112016.txt.gz&dir=data/stdmet/Jan/
%    
%    http://www.ndbc.noaa.gov/download_data.php?filename=4604112016.txt.gz&dir=data/stdmet/Jan/
%    http://www.ndbc.noaa.gov/download_data.php?filename=4604122016.txt.gz&dir=data/stdmet/Feb/
%    
%    sprintf('Reading %s data',num2str(i))
%    [inr,status]=urlread(url);
%    status

%calculate wave energy flux (use equation from wikipedia here)
StdMetData.wavepower=((1024*(9.8)^2)/(64*pi))*(StdMetData.wvht.^2).*StdMetData.dpd;

% basic plotting
figure('color','w')
subplot(5,1,1)
plot(StdMetData.time,StdMetData.pres,'k.')
ylabel('Pressure (mbar)')
datetickzoom('x')

subplot(5,1,2)
plot(StdMetData.time,StdMetData.airtemp,'k.')
ylabel('Air Temp (C)')
datetickzoom('x')

subplot(5,1,3)
plot(StdMetData.time,StdMetData.watertemp,'k.')
ylabel('Water Temp (C)')
datetickzoom('x')

subplot(5,1,4)
plot(StdMetData.time,StdMetData.wavepower,'k.')
ylabel('Wave Power')
datetickzoom('x')

subplot(5,1,5)
plot(StdMetData.time,StdMetData.windspd,'k.')
ylabel('Wind Spd')
datetickzoom('x')

%% store structure
save(sprintf('NDBCStat_%s_StdMetData',num2str(stname)),'data')

% %% store output and calculate trends and patterns
% clc; s=input('Continue with Storing Data Structure and Calculating Trends and Patterns (y/n)?','s')
% if strncmp(s,'y')
%     
%     %store structure
%     save(sprintf('NDBCStat_%s_StdMetData',num2str(stname)),'data')
%     
%     
% end

% % %daily average
% % dy=floor(min(data.dates)):1:ceil(max(data.dates));
% % for i=1:length(dy)
% %     ind=find(floor(data.dates)==dy(i));
% %     if(length(ind)>1) %this calculates a mean if there is only one value in the day
% %         wvs(i,1)=dy(i);
% %         wvs(i,2)=nanmean(data.wavepower(ind));
% %     else
% %         wvs(i,1)=dy(i);
% %         wvs(i,2)=NaN;
% %     end
% % end
% % 
% % figure('color','w')
% % plot(wvs(:,1),wvs(:,2),'k-','linewidth',2)
% % datetick('x','mm/yy')
% % ylabel('Daily Average Wave Power')