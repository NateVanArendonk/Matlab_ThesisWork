% Load in Hsig MAT file 
clearvars
load('HsigMetrics.mat');
[x,y] = deg2utm(lat,lon);
%% Plot
clf
IM = load('E:\Abbas\PS_COSMOS\Thesis_Modeling\GeoTiffs\Puget Sound\SalishSea_UTM.mat');
imagesc(IM.xm,IM.ym,IM.im)
set(gca,'ydir','normal')
cinds = zeros(length(lon),1);
N_c = 30000;
mycolors = parula(N_c);
hold on 
metric = 'yr_max';
switch metric
    case 'yr_max'
        hs = hs_yr;
    case 'max'
        hs = 'max';
end
for ii = 1:length(lon)
    cind = round(hs(ii)*10000);
    %             cind = round(K.hs_smooth(k)*388);
    if cind <= 0
        cind = 1;
    end
    if cind >= 30000
        cind = 30000;
    end
    cinds(ii) = cind;
    plot(x(ii),y(ii), 'o', 'Color', mycolors(cind,:),'MarkerFaceColor',mycolors(cind,:),'MarkerSize',1.8)
    hold on
end

% Adjust color bar
chan = colorbar('Location','EastOutside');
set(chan,'XTick',0:(1/6):1,'XTickLabel',0:.5:3)
chan.FontSize = 18;
chan.FontWeight = 'bold';
switch metric 
    case 'max'
        ylabel(chan,'Maximum Significant Wave Height [m]','FontSize',18,'fontweight','bold')
    case 'yr_max'
        ylabel(chan,'Annual Maximum Significant Wave Height [m]','FontSize',18,'fontweight','bold')
end        
ylabel('Northing [m]','FontSize',18,'fontweight','bold')
xlabel('Easting [m]','FontSize',18,'fontweight','bold')

axis equal
ax = gca;
% Entire Domain 
ax.XLim = [(3.5*10^5),(5.7*10^5)];
ax.YLim = [(5.21*10^6),(5.445*10^6)];

saveFile = 1;
if saveFile
    switch metric
        case 'max'
            printFig(gcf,'MaxHsig_utm',[11 11],'png',300)
        case 'yr_max'
            printFig(gcf,'YEARMaxHsig_utm',[11 11],'png',300)
    end
end