% Load in Hsig MAT file 
clearvars 
clc
H = load('HsigMax.mat');
B = load('HsigMax_Blaine.mat');

%% Grab Unique points so that there is no overlap between the two models 
x = [H.x;B.x];
y = [H.y;B.y];
hs = [H.max_hs;B.max_hs];
temp = [x,y];
[C,ia,ic] = unique(temp,'rows');
hs = hs(ia);
x = C(:,1);
y = C(:,2);

%% Plot
clf
IM = load('E:\Abbas\PS_COSMOS\Thesis_Modeling\GeoTiffs\Whatcom\whatcom1.mat');
imagesc(IM.xm,IM.ym,IM.im)
set(gca,'ydir','normal')
cinds = zeros(length(x),1);
N_c = 30000;
mycolors = parula(N_c);
hold on 
for ii = 1:length(x)
    cind = round(hs(ii)*10000);
    cinds(ii) = cind;
    %             cind = round(K.hs_smooth(k)*388);
    if cind <= 0
        cind = 1;
    end
    if cind >= 30000
        cind = 30000;
    end
    plot(x(ii),y(ii), 'o', 'Color', mycolors(cind,:),'MarkerFaceColor',mycolors(cind,:),'MarkerSize',4)
    hold on
end

% Adjust color bar
chan = colorbar('Location','EastOutside');
set(chan,'XTick',0:(1/3):1,'XTickLabel',0:1:3)
ylabel(chan,'Maximum Significant Wave Height [m]','FontSize',14)
ylabel('Northing [m]','FontSize',14)
xlabel('Easting [m]','FontSize',14)

axis equal
ax = gca;
% Entire Domain 
ax.XLim = [(4.8*10^5),(5.5*10^5)];
ax.YLim = [(5.36*10^6),(5.429*10^6)];

printFig(gcf,'MaxHsig_utm_Whatcom',[11 11],'png',300)