clearvars
addpath C:\Functions_Matlab
addpath C:\Functions_Matlab\mapping\kml

% Load KML with Points
kml = dir('C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\WaveAnalysis\Sensitivity\KML_Points\*.kml');
kmlfol = 'C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\WaveAnalysis\Sensitivity\KML_Points\';
format = [1 4 2 3];

% Load in data to find wave grid point closest to kml in loop below 
fname = 'southern_t200_s10_d300.mat';
D = load(['C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\Hindcast\Southern_LUT_Data\' fname]);

for ii = 1:length(kml)
    t = kml2struct([kmlfol kml(format(ii)).name]);
    K(ii).lon = t.Lon;
    K(ii).lat = t.Lat;
    K(ii).name = t.Name;
    [K(ii).x,K(ii).y] = deg2utm(t.Lat,t.Lon);
    [K(ii).xg,K(ii).yg] = findNearestGridPoint(D.X,D.Y,K(ii).lon,K(ii).lat); % Find closest wave grid point
end

%%
% Specify Dimensions of Plot
p_left = .05;
p_right = .04;
p_top = .03;
p_bot = .07;
p_spacing = .02;
p_wid = ((1-p_right-p_left-p_spacing)/3);
p_height = (1-p_top-p_bot-p_spacing)/3;

% Starting position of plots 
row = 2;
col = 0;

% Speeds to use
speeds = [10 18 25];
ss = 1;
% Tides to use
tides = [-2 1 3];
tt = 1;
clf

% -------------- Direction Sensitivity 
% direction = 'north';
b = parula(5);
dirWant = 0;
tic
for p = 1:9
    axes('position',[p_left+col*(p_wid+p_spacing) p_bot+row*(p_height+p_spacing) p_wid p_height])
    tide2load = 100*tides(tt); % Grab tide
    spd2load = speeds(ss); % Grab speed
    for ii = 1:length(K)
        hsig = zeros(1,length(0:15:345));
        for jj = 1:length(hsig)
            dir2use = wrap2360(90-(dirWant - 180));
            fname = sprintf('southern_t%d_s%d_d%d.mat',tide2load,spd2load,dir2use);
            D = load(['C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\Hindcast\Southern_LUT_Data\' fname]);
            hsig(jj) = D.hs(K(ii).xg,K(ii).yg); % Grab wave point data
            dirWant = dirWant + 15; %Increase direction
        end
        x_axis = 0:15:345;
        ll(ii) = plot(x_axis,hsig,'o','MarkerSize',9,'MarkerFaceColor',b(ii,:),'MarkerEdgeColor',b(ii,:));
        hold on
        line(x_axis,hsig,'Color',b(ii,:));
        hold on
        clear hsig
        dirWant = 0;
    end
    if p == 2 || p == 3 || p == 5 || p == 6
        xlabel(' ')
        ylabel(' ')
        set(gca,'XTickLabel',[])
        set(gca,'YTickLabel',[])
%         xlim([0 18])
%         ylim([0 18])
    elseif p == 1 || p == 4
        xlabel(' ')
        ylabel('Hsig [m]')
        set(gca,'XtickLabel',[])
%         xlim([0 18])
%         ylim([0 18])
    elseif p == 8 || p == 9
        ylabel(' ')
        xlabel('Direction [deg]')
        set(gca,'YtickLabel',[])
%         xlim([0 18])
%         ylim([0 18])
    else
        xlabel('Direction [deg]')
        ylabel('Hsig [m]')
%         xlim([0 18])
%         ylim([0 18])
    end
    % Change to a new row column pair and change the speed and tide level
    col = col + 1;
    ss = ss + 1;
    if col > 2
        col = 0;
        row = row - 1;
        tt = tt + 1; % Use a new tide level 
        ss = 1; % Go back to using the fist speed
    end
    grid on 
    box on
    ylim([0 3.2])
    if p == 1 % add legend 
        lgd = legend(ll,arrayfun(@(x)(x.name),K,'un',0),...
        'interpreter','none',...
        'location','northwest');
        lgd.FontSize = 12;
    end
end
toc
grid on 
box on 

% xlabel('Wind Direction [degrees]','FontSize',12)
% ylabel('Hsig [m]','FontSize',12)
% ax = gca;
% ax.FontSize = 11;
% ax.XLim = [0 350];


