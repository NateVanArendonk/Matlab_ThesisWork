% Explore Swell data from Sean
clearvars
clc

sFol = 'Swell_ModelRuns'; % Location of Swell Runs
d = dir([sFol '\*.mat']); % Get List of all swell runs
load('E:\Abbas\PS_COSMOS\Salish_Model_Resources\WA_Spatial_Data\WA_coast.mat');



MOS = cell(12,2);

D = cell(length(d),12);
tic
for ii = 1:length(d)
    S = load([sFol '\' d(ii).name]); % Load a single run for examination
    
    % Subsample time vector to be just the year that was loaded in
    yearID = strfind(d(ii).name,'_');
    yearID = d(ii).name(yearID(3)+1:yearID(3)+4);
    yearID = str2double(yearID);
    yrInds = year(S.time) == yearID;
    time = S.time(yrInds);
    
    for mm = 1:12
        D{ii,mm} = cell(1,2);
        inds = month(time) == mm;
        D{ii,mm}{1,1} = time(inds);
        D{ii,mm}{1,2} = S.hs_swell(:,inds);
    end
    if ii == length(d)
        lat = S.lat;
        lon = S.lon;
    end
    clear S
    ii
    
end
toc

[x,y] = deg2utm(lat,lon);

%% Calculate statistics about each month for all points

moMn = cell(12,1);
moMx = moMn;
for ii = 1:12
    mons = cell(1,31);
    for yy = 1:length(D)
        mons{yy} = D{yy,ii}{1,2};
    end
    monAvg = cellfun(@(x) mean(x,2),mons,'UniformOutput',false);
    monMax = cellfun(@(x) max(x,[],2),mons,'UniformOutput',false);
    moMn{ii} = mean(cat(3,monAvg{:}),3); % Calculate average monthly mean of this month over all the years
    moMx{ii} = mean(cat(3,monMax{:}),3); % Calculate average monthly max of this month over all the years
end

%% Get Max max and Max mean from lists 
mnMAX = max(cell2mat(moMn));
mxMAX = max(cell2mat(moMx));
return

%% Plot
clf
p_left = .038;
p_right = .05;
p_top = .1;
p_bot = .06;
p_spacing = .01;
p_wid = (1-p_right-p_left-p_spacing)/3;
p_height = (1-p_top-p_bot-p_spacing)/4;

IM = load('E:\Abbas\PS_COSMOS\Thesis_Modeling\GeoTiffs\Puget Sound\SalishSea_UTM.mat');
row = 3;
col = 0;

MONTHS = {'Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'};

metric = 'mean';
tic
for ii = 1:12
    
    axes('position',[p_left+col*(p_wid+p_spacing) p_bot+row*(p_height+p_spacing) p_wid p_height])
    imagesc(IM.xm,IM.ym,IM.im)
    set(gca,'ydir','normal')
    axis equal
    hold on 
    
    col = col + 1;
    if col >= 3
        col = 0;
        row = row - 1;
    end
    
%     axis equal
%     ax = gca;
    xlim([(3.52*10^5),(5.5*10^5)]);
    ylim([(5.3*10^6),(5.4*10^6)]);
    
    
    noXY = [2,3,5,6,8,9];
    noX = [1,4,7];
    noY = [11,12];
    if ismember(ii,noXY)
        set(gca,'XTickLabels',[])
        set(gca,'YTickLabels',[])
    elseif ismember(ii,noX)
        set(gca,'XTickLabels',[])
    elseif ismember(ii,noY)
        set(gca,'YTickLabels',[])
    else
        set(gca,'FontSize',10)
    end
    
    title(MONTHS{ii})
    
    % Plot Wave Metric 
    switch metric 
        case 'max'
            
            % Colors
            N_c = round(mxMAX)*10000;
            mycolors = parula(N_c);
            
            for mm = 1:length(lon)
                cind = round(moMx{ii,1}(mm)*10000);
                if cind <= 0
                    cind = 1;
                end
                if cind >= N_c
                    cind = N_c;
                end
                plot(x(mm),y(mm), 'o', 'Color', mycolors(cind,:),'MarkerFaceColor',mycolors(cind,:),'MarkerSize',2)
                hold on
            end
            saveNM = 'AvgMonthlyMax';
            
            chan = colorbar('Location','NorthOutside');
            set(chan,'Position',[p_left .93 ((p_left + 3*p_wid + 2*p_spacing)-p_left) .02])
            set(chan,'XTick',0:(1/4):1,'XTickLabel',0:1:4)
            chan.Label.String = 'Hs [m]';
            chan.FontSize = 10;
            
        case 'mean'
            % Colors
            N_c = round(mnMAX)*10000;
            mycolors = parula(N_c);
            
            for mm = 1:length(lon)
                cind = round(moMn{ii,1}(mm)*10000);
                if cind <= 0
                    cind = 1;
                end
                if cind >= N_c
                    cind = N_c;
                end
                plot(x(mm),y(mm), 'o', 'Color', mycolors(cind,:),'MarkerFaceColor',mycolors(cind,:),'MarkerSize',2)
                hold on
            end
            saveNM = 'AvgMonthlyMean';
            chan = colorbar('Location','NorthOutside');
            set(chan,'Position',[p_left .93 ((p_left + 3*p_wid + 2*p_spacing)-p_left) .02])
            set(chan,'XTick',0:(1/4):1,'XTickLabel',0:.5:2)
            chan.Label.String = 'Hs [m]';
            chan.FontSize = 10;
    end
end
toc


printFig(gcf,saveNM,[11 9],'png',300)
