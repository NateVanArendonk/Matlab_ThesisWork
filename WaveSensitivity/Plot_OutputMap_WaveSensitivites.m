clearvars

addpath C:\Functions_Matlab

fol_run = 'C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\Hindcast\Southern_LUT_Data\';
IM = load('C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\GeoTiffs\South_Central_Entire.mat');

dir = 'north';
spd = 10;

% Specify Dimensions of Plot
p_left = .05;
p_right = .05;
p_top = .1;
p_bot = .06;
p_spacing = .02;
p_wid = ((1-p_right-p_left-p_spacing)/3);
p_height = ((1-p_top-p_bot-p_spacing)/2)+.02;

row = 1; 
col = 0;
clf
tic
% Plot speed sensitivity
for ii = 1:6
    % load data
    spd2load = round(spd);
%     fprintf('Should be loading %d',spd);
%     fprintf('  - Loading Speed: %d\n',spd2load);
    switch dir 
        case 'north'
            fname = sprintf('southern_t200_s%d_d300.mat',spd2load);
        case 'south'
            fname = sprintf('southern_t200_s%d_d60.mat',spd2load);
    end
    data = load([fol_run fname]);
    [X_utm,Y_utm,~] = deg2utm(data.Y(:),data.X(:));
    [r,c] = size(data.X);
    X_utm = reshape(X_utm,r,c); Y_utm = reshape(Y_utm,[r,c]);
    inds = find(data.depth <= 0.5);
    data.hs(inds) = NaN;
    spd = spd+2.5;
    
    
    % Plot
    axes('position',[p_left+col*(p_wid+p_spacing) p_bot+row*(p_height+p_spacing) p_wid p_height])
    imagesc(IM.xm,IM.ym,IM.im)
    set(gca,'ydir','normal')
    hold on
    pcolor(X_utm,Y_utm,data.hs)
    shading interp
%     colorbar
    caxis([0 2.5])
    axis equal
%     chan = colorbar('Location','EastOutside');
%     ylabel(chan,'Significant Wave Height [m]')  
    ax = gca;
    
    ax.XLim = [(4.88*10^5), (5.62*10^5)];
    ax.YLim = [(5.211*10^6), (5.34*10^6)];
    
    
    col = col + 1;
    if col > 2
        col = 0;
        row = row - 1;
    end
    if ii == 1 
        set(gca,'XTickLabel',[])
        ylabel('Northing [m]')
        xlabel(' ')
    elseif ii == 2 || ii == 3 
        set(gca,'XTickLabel',[])
        set(gca,'YTickLabel',[])
        ylabel(' ')
        xlabel(' ')
    elseif ii == 5 || ii == 6
        set(gca,'YTickLabel',[])
        ylabel(' ')
        xlabel('Easting [m]')
    else
        xlabel('Easting [m]')
        ylabel('Northing [m]')
    end
    
end
chan = colorbar('Location','NorthOutSide','FontSize',8);
set(chan,'Position',[.03 1-.04 ((p_left + 3*p_wid + 3*p_spacing)-p_left) .02])
toc
%%
printFig(gcf,'ChangingWindSpeedSnesitivity',[11 8.5],'png',300)
%% Plot Directional Sensitivity 

% Specify Dimensions of Plot
p_left = .05;
p_right = .05;
p_top = .1;
p_bot = .06;
p_spacing = .02;
p_wid = ((1-p_right-p_left-p_spacing)/4);
p_height = ((1-p_top-p_bot-p_spacing)/2)+.02;

row = 1; 
col = 0;
clf
tic
dirWant = 120;
direction = wrap2360(90-(dirWant - 180));
for ii = 1:8
    % load data
    spd2load = 23;

    fname = sprintf('southern_t200_s%d_d%d.mat',spd2load,direction);
    
    data = load([fol_run fname]);
    [X_utm,Y_utm,~] = deg2utm(data.Y(:),data.X(:));
    [r,c] = size(data.X);
    X_utm = reshape(X_utm,r,c); Y_utm = reshape(Y_utm,[r,c]);
    inds = find(data.depth <= 0.5);
    data.hs(inds) = NaN;
    direction = direction + 15;

    % Plot
    axes('position',[p_left+col*(p_wid+p_spacing) p_bot+row*(p_height+p_spacing) p_wid p_height])
    imagesc(IM.xm,IM.ym,IM.im)
    set(gca,'ydir','normal')
    hold on
    pcolor(X_utm,Y_utm,data.hs)
    shading interp
%     colorbar
    caxis([0 2.5])
    axis equal
%     chan = colorbar('Location','EastOutside');
%     ylabel(chan,'Significant Wave Height [m]')  
    ax = gca;
    
    ax.XLim = [(4.88*10^5), (5.62*10^5)];
    ax.YLim = [(5.211*10^6), (5.34*10^6)];
    
    
    col = col + 1;
    if col > 3
        col = 0;
        row = row - 1;
    end
    if ii == 1 
        set(gca,'XTickLabel',[])
        ylabel('Northing [m]')
        xlabel(' ')
    elseif ii == 2 || ii == 3 || ii == 4
        set(gca,'XTickLabel',[])
        set(gca,'YTickLabel',[])
        ylabel(' ')
        xlabel(' ')
    elseif ii == 6 || ii == 7 || ii == 8
        set(gca,'YTickLabel',[])
        ylabel(' ')
        xlabel('Easting [m]')
    else
        xlabel('Easting [m]')
        ylabel('Northing [m]')
    end
    
end
chan = colorbar('Location','NorthOutSide','FontSize',8);
set(chan,'Position',[.069 1-.04 ((p_left + 4*p_wid + 3*p_spacing)-p_left-.069) .02])
toc

%%
printFig(gcf,'ChangingWindDirectionSnesitivity',[11 8.5],'png',300)
