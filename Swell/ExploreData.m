% Explore Swell data from Sean
clearvars 
clc

sFol = 'Swell_ModelRuns'; % Location of Swell Runs
d = dir([sFol '\*.mat']); % Get List of all swell runs
avgSwell = NaN(6352,31);
% Output is a 3 hourly time vector 
return
tic % 6 minutes 
for ii = 1:length(d)
    S = load([sFol '\' d(ii).name]); % Load a single run for examination
%     yearID = strfind(d(ii).name,'_');
%     yearID = d(ii).name(yearID(3)+1:yearID(3)+4);
%     yearID = str2double(yearID);
%     yrInds = year(S.time) == yearID;
%     time = S.time(yrInds);
    avgSwell(:,ii) = max(S.hs_swell,[],2);
    clear S 
end
toc

yrAvgMax = mean(avgSwell,2);
% Note time variable contains time from the entire run not for that single
% year 
[x,y] = deg2utm(S.lat,S.lon);
%% Plot 
clf
IM = load('E:\Abbas\PS_COSMOS\Thesis_Modeling\GeoTiffs\Puget Sound\SalishSea_UTM.mat');
imagesc(IM.xm,IM.ym,IM.im)
set(gca,'ydir','normal')
cinds = zeros(length(x),1);
N_c = 40000;
mycolors = parula(N_c);
hold on 

for ii = 1:length(S.lat)
    cind = round(yrAvgMax(ii)*10000);
    cinds(ii) = cind;
    if cind <= 0
        cind = 1;
    end
    if cind >= N_c
        cind = N_c;
    end
    plot(x(ii),y(ii), 'o', 'Color', mycolors(cind,:),'MarkerFaceColor',mycolors(cind,:),'MarkerSize',4)
    hold on
    
end

% Adjust color bar
chan = colorbar('Location','EastOutside');
set(chan,'XTick',0:(1/4):1,'XTickLabel',0:1:4)
ylabel(chan,'Average Annual Max Swell Wave Height [m]','FontSize',14)
ylabel('Northing [m]','FontSize',14)
xlabel('Easting [m]','FontSize',14)

axis equal
ax = gca;
ax.XLim = [(3.5*10^5),(5.7*10^5)];
ax.YLim = [(5.21*10^6),(5.4279*10^6)];
%% 
printFig(gcf,'SwellAvgAnnualMax',[11 11],'png',300)