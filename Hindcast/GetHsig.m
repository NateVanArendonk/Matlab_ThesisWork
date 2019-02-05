% From All Hindcast Timeseries files, Grab Hsig
clearvars

% Location and list of all hindcast files
hindcastFolder = 'E:\Abbas\PS_COSMOS\Thesis_Modeling\LUT\Hindcast_Output\';
h = dir('E:\Abbas\PS_COSMOS\Thesis_Modeling\LUT\Hindcast_Output\*.mat');

% Get names of all the stations - Just station name, no unique ID
names = strings(length(h),1);
for ii = 1:length(h)
    [~,tempName,~] = fileparts(h(ii).name); % grab name from file
    underInds = strfind(tempName,'_'); % Find where underscores are
    names(ii) = tempName(1:underInds(1)-1); % Get rid of underscores from name
end
uniqueNames = unique(names); % Get a list of each unique station
%% Go in and load each part of the basin and get statistics of choice
tic
max_hs = []; % Max Hsig from Time Series
lat = []; % Latitude
lon = []; % Longitude
x = []; % X Utm
y = []; % Y Utm
hs_yr = [];

% Load a run to get time vector to grab yearly indices 
temp = load([hindcastFolder h(1).name]); % smallest file 
time = temp.time;
twl = temp.twl;
years = year(time(1)):1:year(time(end));
yrInds = cell(length(years),1);
for ii = 1:length(yrInds)
    yrInds{ii} = find(year(time) == years(ii));
end

% To house conditions that bring max waves
maxHstwl = [];
maxHsspeed = [];
maxHswnddir = [];
maxHstime = [];
maxHsdir = [];
for ii = 1:length(uniqueNames)
    basinName = uniqueNames(ii); % Name of basin
    basinInds = strcmp(basinName,names); % Find which indices are correct basin
    basinInds = find(basinInds == 1); % Get exact indice
    for jj = 1:length(basinInds) % now load each part of basin
        B = load([hindcastFolder h(basinInds(jj)).name]);
        % Calculate Hsig from time series
        [tempMax,I] = max(B.hs_ts,[],2);
        max_hs = vertcat(max_hs,tempMax);
        maxHsspeed = vertcat(maxHsspeed,B.speed(I));
        maxHswnddir = vertcat(maxHswnddir,B.wnddir(I));
        maxHstime = vertcat(maxHstime,B.time(I)');
        maxHstwl = vertcat(maxHstwl,B.twl(I));
        maxHsdir = vertcat(maxHsdir,B.hdir(I));
        lon = vertcat(lon,B.lon);
        lat = vertcat(lat,B.lat);
%         x = vertcat(x,B.x_u);
%         y = vertcat(y,B.y_u);
        % Calculate yearly max 
        temp_yr = zeros(length(B.lon),length(years));
        for yy = 1:length(yrInds)
            inds = yrInds{yy};
            hs = B.hs_ts;
            hs_yr_temp = hs(:,inds);
            [temp_yr(:,yy), ~] = max(hs_yr_temp,[],2);
        end
        temp_yr = mean(temp_yr,2);
        hs_yr = vertcat(hs_yr,temp_yr);
    end
    fprintf('Completed %s Basin...Moving On\n',basinName);
    fprintf('%2.1f Percent Complete\n',ii/length(uniqueNames)*100);
    
end
toc
save('HsigMetrics','lat','lon','max_hs','hs_yr','maxHsspeed','time','maxHstime','maxHstwl','twl','maxHswnddir','maxHsdir')
return
%% Plot
clf
IM = load('E:\Abbas\PS_COSMOS\Thesis_Modeling\GeoTiffs\Puget Sound\SalishSea_UTM.mat');
imagesc(IM.xm,IM.ym,IM.im)
set(gca,'ydir','normal')
cinds = zeros(length(lon),1);
N_c = 50000;
mycolors = parula(N_c);
hold on 
for ii = 1:length(lon)
    cind = round(max_hs(ii)*10000);
    cinds(ii) = cind;
    %             cind = round(K.hs_smooth(k)*388);
    if cind <= 0
        cind = 1;
    end
    plot(x(ii),y(ii), 'o', 'Color', mycolors(cind,:),'MarkerFaceColor',mycolors(cind,:),'MarkerSize',4)
    hold on
end

% Adjust color bar
chan = colorbar('Location','EastOutside');
set(chan,'XTick',0:(1/5):1,'XTickLabel',0:1:5)
ylabel(chan,'Maximum Significant Wave Height [m]','FontSize',14)
ylabel('Northing [m]','FontSize',14)
xlabel('Easting [m]','FontSize',14)

axis equal
ax = gca;
ax.XLim = [(3.5*10^5),(5.7*10^5)];
ax.YLim = [(5.21*10^6),(5.4279*10^6)];
% printFig(gcf,'MaxHsig_utm',[11 11],'png',300)