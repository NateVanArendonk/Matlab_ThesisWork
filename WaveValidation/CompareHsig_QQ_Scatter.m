% Load in wave hindcast output
clearvars

w = dir('Hindcast_Output/*.mat');
nameVec = strings(length(w),1);
for ii = 1:length(w)
    W(ii) = load(['Hindcast_Output/' w(ii).name]);
    nameVec(ii) = strcat('H',num2str(ii));
end
H = load('E:\Abbas\Model_Met_Forcings\HeinBank\obs_dungeness2006.mat');
clear h ii w


% Get obs for same time period as model
% Note that the only good agreement time looks to be in 2004 so I'll
% subsample the model output to be 2004 as well as the obs

% Now subsample the model to just be in the window of the obs
modelStart = find(W(1).time == H.time(1));
modelEnd = find(W(1).time == H.time(end));
modelInds = modelStart:modelEnd;
for ii = 1:length(W)
    W(ii).hs_ts = W(ii).hs_ts(modelInds);
    W(ii).speed = W(ii).speed(modelInds);
    W(ii).time = W(ii).time(modelInds);
    W(ii).tm_ts = W(ii).tm_ts(modelInds);
    W(ii).tp_ts = W(ii).tp_ts(modelInds);
    W(ii).twl = W(ii).twl(modelInds);
    W(ii).wnddir = W(ii).wnddir(modelInds);
end

% Interp the model predictions on to the time vector of obs
tvec = H.time;
for ii = 1:length(W)
    W(ii).hs_ts = interp1(W(ii).time,W(ii).hs_ts,tvec);
    W(ii).speed = interp1(W(ii).time,W(ii).speed ,tvec);
    W(ii).tm_ts = interp1(W(ii).time,W(ii).tm_ts,tvec);
    W(ii).tp_ts = interp1(W(ii).time,W(ii).tp_ts,tvec);
    W(ii).twl = interp1(W(ii).time,W(ii).twl,tvec);
    W(ii).wnddir = interp1(W(ii).time,W(ii).wnddir,tvec);
    W(ii).time = tvec;
end

%% MultiPannel QQ Plots 

p_left = .05;
p_right = .05;
p_top = .03;
p_bot = .07;
p_spacing = .01;
p_wid = ((1-p_right-p_left-p_spacing)/3);
p_height = (1-p_top-p_bot-p_spacing)/3;

row = 2;
col = 0;

% hist = find(year(K(1).time) >= 1950 & year(K(1).time) <= 2020);
% fut = find(year(K(1).time) >= 2020 & year(K(1).time) <= 2100);
% hist(1:50) = [];

scenario = '85Maxs';
clf
for ii = 1:length(W)

    cdf1 = sort(H.Hs,'ascend');
    cdf1y = linspace(0,1,length(cdf1));
    
    cdf2 = sort(W(ii).hs_ts,'ascend');
    cdf2y = linspace(0,1,length(cdf2));
    
    %Interp cdfs
    cdf2i = interp1(cdf2y,cdf2,cdf1y);
    
    x_axis = linspace(0,1,length(cdf1));
    axes('position',[p_left+col*(p_wid+p_spacing) p_bot+row*(p_height+p_spacing) p_wid p_height])
    scatter(H.Hs,W(ii).hs_ts)
    hold on
    plot(cdf1,cdf2i,'r+')
    tit = sprintf('%s',nameVec(ii));
    title(tit)
    set(get(gca,'title'),'Position',[3.75 0.36 1.00011])
    grid on
    hold on
    line([0 50],[0 50],'Color','k','LineStyle','--')
    if ii == 2 || ii == 3 || ii == 5 || ii == 6
        xlabel(' ')
        ylabel(' ')
        set(gca,'XTickLabel',[])
        set(gca,'YTickLabel',[])
        xlim([0 4])
        ylim([0 4])
    elseif ii == 1 || ii == 4
        set(gca,'XtickLabel',[])
        xlim([0 4])
        ylim([0 4])
    elseif ii == 8 || ii == 9
        set(gca,'YtickLabel',[])
        xlim([0 4])
        ylim([0 4])
    else  
        xlim([0 4])
        ylim([0 4])
    end
    set(gca,'FontSize',14)
    
    % Increment position 
    col = col + 1;
    if col >=3
        col = 0;
        row = row - 1;
    end
end

printFig(gcf,'HeinVSLUT_QQScatter_Mulitpanel',[14 14],'png',300)