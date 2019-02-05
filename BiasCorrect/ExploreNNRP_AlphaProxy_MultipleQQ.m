clearvars
np = dir('NNRP_Data\*.mat');
np = nestedSortStruct(np,{'name'});
for ii = 1:length(np)
    Nstring = strfind(np(ii).name,'N');
    Pstring = strfind(np(ii).name,'.');
    num = np(ii).name(Nstring+1:Pstring-1);
    np(ii).num = str2num(num);
end
np = nestedSortStruct(np,{'num'});

for ii = 1:length(np)
    N(ii) = load(['NNRP_Data\' np(ii).name]);
    %     N(ii).name = strcat('N',num2str(ii));
end
for ii = 1:length(np)
    N(ii).name = strcat('N',num2str(ii));
end

waterID = 5;

%% Grab Water point and box around it (1 - 9 points)



clf
p_left = .05;
p_right = .05;
p_top = .07;
p_bot = .07;
p_spacing = .01;
p_wid = ((1-p_right-p_left-p_spacing)/3);
p_height = (1-p_top-p_bot-p_spacing)/6;

row = 5;
col = 0;
b = parula(9);

shape = {'+','o','*','x','s','d','^','p','>'};

for ii = 1:length(T)
    cdf1 = sort(T(waterID).wndspd,'ascend');
    cdf1y = linspace(0,1,length(cdf1));
    
    cdf2 = sort(T(waterID).wndspd,'ascend');
    cdf2y = linspace(0,1,length(cdf2));
    
    %Interp cdfs
    cdf2i = interp1(cdf2y,cdf2,cdf1y);
    
    x_axis = linspace(0,1,length(cdf1));
%     scatter(T(5).wndspd,T(ii).wndspd)
    hold on
    set(gca,'Color',[.7 .7 .7])
    cc = randi([1 100],1);
    cc = ii*cc;
    ll(ii) = plot(cdf1,cdf2i,'Color',b(ii,:),'LineWidth',2);
    ll(ii) = plot(cdf1(end),cdf2i(end),shape{ii},'Color',b(ii,:),'MarkerFaceColor',b(ii,:),'MarkerSize',10);

%     text(1,15,tit,'FontSize',14)
    grid on
    hold on
    line([0 50],[0 50],'Color','k','LineStyle','--')
    set(gca,'FontSize',14)
    xticks([0:5:20])
    yticks([0:5:20])
    xlim([0 20])
    ylim([0 20])

end




lgd = legend(ll,arrayfun(@(x)(x.name),N,'un',0),...
    'interpreter','none',...
    'location','northwest');
lgd.FontSize = 12;

printFig(gcf,'QQ_SurroundingBlake_AllDirections',[11 8.5],'png',300)
%% Grab Water point and box around it (1 - 9 points) for different wind regimes


b = parula(length(N));

shape = {'+','o','*','x','s','d','^','p','>'};



row = 5;
col = 0;
wind_direction = 'South';
clf
for ii = 1:length(N)

    % -------------- Base --------------------------
    base = N(waterID).wnddir;
    Ninds = base >= 290 | base <= 45;
    Sinds = base >= 135 & base <= 225;
    switch wind_direction
        case 'North'
            inds = Ninds;
        case 'South'
            inds = Sinds;
    end
    baseSpd = N(waterID).wndspd;
    cdf1 = sort(baseSpd(inds),'ascend');
    cdf1y = linspace(0,1,length(cdf1));
    
    % --------------- Comparison Data -------------------
    compare = N(ii).wnddir;
    Ninds = compare >= 290 | compare <= 45;
    Sinds = compare >= 135 & compare <= 225;
    switch wind_direction
        case 'North'
            inds = Ninds;
        case 'South'
            inds = Sinds;
    end
    compSpd = N(ii).wndspd;
    cdf2 = sort(compSpd(inds),'ascend');
    cdf2y = linspace(0,1,length(cdf2));
    
    %Interp cdfs
    cdf2i = interp1(cdf2y,cdf2,cdf1y);
    
    x_axis = linspace(0,1,length(cdf1));
%     scatter(T(5).wndspd,T(ii).wndspd)
    hold on
    set(gca,'Color',[.7 .7 .7])
    cc = randi([1 100],1);
    cc = ii*cc;
    ll(ii) = plot(cdf1,cdf2i,'Color',b(ii,:),'LineWidth',2);
    ll(ii) = plot(cdf1(end),cdf2i(end),shape{ii},'Color',b(ii,:),'MarkerFaceColor',b(ii,:),'MarkerSize',12);
        line([0 50],[0 50],'Color','k','LineStyle','--')
    grid on 
    box on 
    switch wind_direction
        case 'North'
            xlim([0 13])
            ylim([0 13])
        case 'South'
            xlim([0 18])
            ylim([0 18])
    end
end

lgd = legend(ll,arrayfun(@(x)(x.name),N,'un',0),...
    'interpreter','none',...
    'location','northwest');
lgd.FontSize = 12;
set(gcf,'InvertHardCopy','off')
switch wind_direction
    case 'North'
        printFig(gcf,'QQ_SurroundingBlake_NorthBlows',[15 11],'png',300)
    case 'South'
        printFig(gcf,'QQ_SurroundingBlake_SouthBlows',[15 11],'png',300)
end


%% Make Plot of extracted points with qq plots 
clf
load('C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\Quantile_Correction\WA_coast.mat');
% Find all shapes above threshold
wa_lat = [];
wa_lon = [];
thresh = 500 ;
for j = 1:length(wa_coast)
    temp_x = wa_coast(j).X;
    temp_y = wa_coast(j).Y;
    if length(temp_x) >= thresh %&& j ~= 3348 % 3348 is oregon
        for m = 1:length(temp_x);
            wa_lat(end+1) = temp_y(m);
            wa_lon(end+1) = temp_x(m);
        end
    end
end
b = parula(length(N));
subplot(1,3,1)
plot(wa_lon,wa_lat)
hold on 
shape = {'+','o','*','x','s','d','^','p','>'};
for ii = 1:length(N)
    pp(ii) = plot(N(ii).lon,N(ii).lat,shape{ii},'Color',b(ii,:),'MarkerFaceColor',b(ii,:),'MarkerSize',12);
    hold on 
end

xlim([-122.75 -122.25])
ylim([47.1 47.9])


lgd = legend(pp,arrayfun(@(x)(x.name),N,'un',0),...
    'interpreter','none',...
    'location','northwest');
lgd.FontSize = 12;


subplot(1,3,2:3)


wind_direction = 'South';
for ii = 1:length(N)

    % -------------- Base --------------------------
    base = N(waterID).wnddir;
    Ninds = base >= 290 | base <= 45;
    Sinds = base >= 135 & base <= 225;
    switch wind_direction
        case 'North'
            inds = Ninds;
        case 'South'
            inds = Sinds;
    end
    baseSpd = N(waterID).wndspd;
    cdf1 = sort(baseSpd(inds),'ascend');
    cdf1y = linspace(0,1,length(cdf1));
    
    % --------------- Comparison Data -------------------
    compare = N(ii).wnddir;
    Ninds = compare >= 290 | compare <= 45;
    Sinds = compare >= 135 & compare <= 225;
    switch wind_direction
        case 'North'
            inds = Ninds;
        case 'South'
            inds = Sinds;
    end
    compSpd = N(ii).wndspd;
    cdf2 = sort(compSpd(inds),'ascend');
    cdf2y = linspace(0,1,length(cdf2));
    
    %Interp cdfs
    cdf2i = interp1(cdf2y,cdf2,cdf1y);
    
    x_axis = linspace(0,1,length(cdf1));
%     scatter(T(5).wndspd,T(ii).wndspd)
    hold on
    set(gca,'Color',[.7 .7 .7])
    cc = randi([1 100],1);
    cc = ii*cc;
    ll(ii) = plot(cdf1,cdf2i,'Color',b(ii,:),'LineWidth',2);
    ll(ii) = plot(cdf1(end),cdf2i(end),shape{ii},'Color',b(ii,:),'MarkerFaceColor',b(ii,:),'MarkerSize',12);
        line([0 50],[0 50],'Color','k','LineStyle','--')
    grid on 
    box on 
    switch wind_direction
        case 'North'
            xlim([0 13])
            ylim([0 13])
        case 'South'
            xlim([0 18])
            ylim([0 18])
    end
    axis equal
end

lgd = legend(ll,arrayfun(@(x)(x.name),N,'un',0),...
    'interpreter','none',...
    'location','northwest');
lgd.FontSize = 12;
set(gcf,'InvertHardCopy','off')
switch wind_direction
    case 'North'
        printFig(gcf,'QQ_SurroundingBlake_NorthBlows_Spatial',[15 11],'png',300)
    case 'South'
        printFig(gcf,'QQ_SurroundingBlake_SouthBlows_Spatial',[15 11],'png',300)
end
%% Plot North and South 
clf
load('C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\Quantile_Correction\WA_coast.mat');
% Find all shapes above threshold
wa_lat = [];
wa_lon = [];
thresh = 500 ;
for j = 1:length(wa_coast)
    temp_x = wa_coast(j).X;
    temp_y = wa_coast(j).Y;
    if length(temp_x) >= thresh %&& j ~= 3348 % 3348 is oregon
        for m = 1:length(temp_x);
            wa_lat(end+1) = temp_y(m);
            wa_lon(end+1) = temp_x(m);
        end
    end
end
b = parula(length(N));
subplot(1,3,1)
plot(wa_lon,wa_lat)
hold on 
shape = {'+','o','*','x','s','d','^','p','>'};
for ii = 1:length(N)
    pp(ii) = plot(N(ii).lon,N(ii).lat,shape{ii},'Color',b(ii,:),'MarkerFaceColor',b(ii,:),'MarkerSize',12);
    hold on 
end

xlim([-122.85 -122.45])
ylim([47 47.5])
lgd = legend(pp,arrayfun(@(x)(x.name),N,'un',0),...
    'interpreter','none',...
    'location','northwest');
lgd.FontSize = 12;



subplot(1,3,2:3)
for ii = 1:length(N)

    % -------------- Base --------------------------
    % ----------- North First ----------------------
    base = N(waterID).wnddir;
    Ninds = base >= 290 | base <= 45;
    Sinds = base >= 135 & base <= 225;
    baseSpd = N(waterID).wndspd;
    cdf1 = sort(baseSpd(Ninds),'ascend');
    cdf1y = linspace(0,1,length(cdf1));
    
    % --------------- Comparison Data -------------------
    compare = N(ii).wnddir;
    Ninds = compare >= 290 | compare <= 45;
    Sinds = compare >= 135 & compare <= 225;
    compSpd = N(ii).wndspd;
    cdf2 = sort(compSpd(Ninds),'ascend');
    cdf2y = linspace(0,1,length(cdf2));
    %Interp cdfs
    cdf2i = interp1(cdf2y,cdf2,cdf1y);
    x_axis = linspace(0,1,length(cdf1));

    hold on
    set(gca,'Color',[.7 .7 .7])
    plot(cdf1,cdf2i,'Color',b(ii,:),'LineWidth',2);
    nn(ii) = plot(cdf1(end),cdf2i(end),shape{ii},'Color',b(ii,:),'MarkerFaceColor',b(ii,:),'MarkerSize',12);
        line([0 50],[0 50],'Color','k','LineStyle','--')

        
        % -------------- Base --------------------------
    % ----------- North First ----------------------
    base = N(waterID).wnddir;
    Ninds = base >= 290 | base <= 45;
    Sinds = base >= 135 & base <= 225;
    baseSpd = N(waterID).wndspd;
    cdf1 = sort(baseSpd(Sinds),'ascend');
    cdf1y = linspace(0,1,length(cdf1));
    
    % --------------- Comparison Data -------------------
    compare = N(ii).wnddir;
    Ninds = compare >= 290 | compare <= 45;
    Sinds = compare >= 135 & compare <= 225;
    compSpd = N(ii).wndspd;
    cdf2 = sort(compSpd(Sinds),'ascend');
    cdf2y = linspace(0,1,length(cdf2));
    %Interp cdfs
    cdf2i = interp1(cdf2y,cdf2,cdf1y);
    x_axis = linspace(0,1,length(cdf1));

    hold on
    set(gca,'Color',[.7 .7 .7])
    plot(cdf1,cdf2i,'Color',b(ii,:),'LineWidth',2,'LineStyle','--');
    ss(ii) = plot(cdf1(end),cdf2i(end),shape{ii},'Color',b(ii,:),'MarkerFaceColor',b(ii,:),'MarkerSize',12);
    line([0 50],[0 50],'Color','k','LineStyle','--')

    grid on 
    box on 
    axis equal
end
xlim([0 18])
ylim([0 18])
xticks([0:2:18])
yticks([0:2:18])

lgdN = legend(nn,arrayfun(@(x)(x.name),N,'un',0),...
    'interpreter','none',...
    'location','northwest');
lgdN.FontSize = 12;

lgdS = legend(ss,arrayfun(@(x)(x.name),N,'un',0),...
    'interpreter','none',...
    'location','northwest');
lgdS.FontSize = 12;

set(gcf,'InvertHardCopy','off')
printFig(gcf,'QQ_Spatial_NorthAndSouth_MultipleStations_SouthTacoma',[15 8.5],'png',300)


