% Load in wave hindcast output
clearvars

M = load('Hindcast_Output/H5_waveout.mat');

% Load in Buoy Data
% Load in Hein Bank Data
h = dir('E:\Abbas\Model_Met_Forcings\HeinBank\*.mat');
inds2get = 2:5;
H = struct;
H.time = [];
H.Hs = [];
H.Tp = [];

% Load 2005 - 2008 waves 
for ii = 1:length(inds2get) % 2005 - 2008
    temp = load(['E:\Abbas\Model_Met_Forcings\HeinBank\' h(inds2get(ii)).name]);
    H.time = vertcat(H.time,temp.time);
    H.Hs = vertcat(H.Hs,temp.Hs);
    H.Tp = vertcat(H.Tp,temp.peakT);
end
clear h ii w inds2get
O = H;
clear H;

time = datenum(2005,1,1,00,00,00):1/24:datenum(2008,12,31,23,00,00);
O.hs = interp1(O.time,O.Hs,time);
O.time_i = time;
return
%%
%Get obs for same time period as model

% Note that the only good agreement time looks to be in 2004 so I'll
% subsample the model output to be 2004 as well as the obs

% Now subsample the model to just be in the window of the obs
modelStart = find(M.time == O.time_i(1));
modelEnd = find(M.time == O.time_i(end));
modelInds = modelStart:modelEnd;
for ii = 1:length(M)
    M(ii).hs_ts = M(ii).hs_ts(modelInds);
    M(ii).speed = M(ii).speed(modelInds);
    M(ii).time = M(ii).time(modelInds);
    M(ii).tm_ts = M(ii).tm_ts(modelInds);
    M(ii).tp_ts = M(ii).tp_ts(modelInds);
    M(ii).twl = M(ii).twl(modelInds);
    M(ii).wnddir = M(ii).wnddir(modelInds);
end



% Find NaN in obs and make model NaN
nInds = isnan(O.hs);
M.hs_ts(nInds) = NaN;
M.tp_ts(nInds) = NaN;
% Interp the model predictions on to the time vector of obs
% tvec = D.time_i;
% for ii = 1:length(W)
%     W(ii).hs_ts = interp1(W(ii).time,W(ii).hs_ts,tvec);
%     W(ii).speed = interp1(W(ii).time,W(ii).speed ,tvec);
%     W(ii).tm_ts = interp1(W(ii).time,W(ii).tm_ts,tvec);
%     W(ii).tp_ts = interp1(W(ii).time,W(ii).tp_ts,tvec);
%     W(ii).twl = interp1(W(ii).time,W(ii).twl,tvec);
%     W(ii).wnddir = interp1(W(ii).time,W(ii).wnddir,tvec);
%     W(ii).time = tvec;
% end


%% Plot Timeseries with RMSE and Bias 
%
clf
o = plot(O.time_i,O.hs);
datetick()
hold on 
p = plot(M.time,M.hs_ts,'o','color',[.7 .7 .7]);
lgd = legend([o,p],'Hein Bank Obs','LUT Predictions','Location','NorthEast');
lgd.FontSize = 12;
ylabel('Wave Height [m]')
xlabel('Time [years]')
set(gca,'FontSize',10)
grid on
rmse = sqrt(nanmean((M.hs_ts - O.hs').^2));
rmse = sprintf('RMSE: %.2f',rmse);
text(datenum(2005,5,1),4.0,rmse,'FontSize',15)
bias = nanmean(M(1).hs_ts - O.hs');
bias = sprintf('Bias: %.2f',bias);
text(datenum(2005,5,1),3.9,bias,'FontSize',15);
% xticks([datenum(2005,1,1) datenum(2006,1,1) datenum(2007,1,1) datenum(2008,1,1)]);
% xticklabels({'2005','2006','2007','2008'})
set(gca,'FontSize',10)
printFig(gcf,'HeinBank_NNRP_RMSE',[11 8.5],'png')
%% Compare binned wind speeds and how well wave heights are caputred in each bin -- NEEDS EDITING 
b1 = M.speed > 0 & M.speed <= 5;
b2 = M.speed > 5 & M.speed <= 10;
b3 = M.speed > 10 & M.speed <= 15;
b4 = M.speed > 15;
clf
% 0 - 5 m/s bin
subplot(4,1,1)
o = plot(O.time_i(b1),O.hs(b1));
datetick()
hold on 
p = plot(M.time(b1),M.hs_ts(b1),'o','color',[.7 .7 .7]);
lgd = legend([o,p],'Sentry Shoal Obs','LUT Predictions','Location','NorthEast');
lgd.FontSize = 12;
ylabel('Wave Height [m]')
set(gca,'FontSize',10)
grid on
rmse = sqrt(nanmean((M.hs_ts(b1) - O.hs(b1)).^2));
rmse = sprintf('RMSE: %.2f',rmse);
text(datenum(1993,5,1),2.5,rmse,'FontSize',15)
bias = nanmean(M(1).hs_ts(b1) - O.hs(b1));
bias = sprintf('Bias: %.2f',bias);
text(datenum(1993,5,1),2.1,bias,'FontSize',15);
title('Wind speeds: 0 - 5 m/s')

% 5 - 10 m/s bin
subplot(4,1,2)
o = plot(D.time_i(b2),D.hs(b2));
datetick()
hold on 
p = plot(W.time(b2),W.hs_ts(b2),'o','color',[.7 .7 .7]);
% lgd = legend([o,p],'Sentry Shoal Obs','LUT Predictions','Location','NorthEast');
% lgd.FontSize = 12;
ylabel('Wave Height [m]')
set(gca,'FontSize',10)
grid on
rmse = sqrt(nanmean((W.hs_ts(b2) - D.hs(b2)).^2));
rmse = sprintf('RMSE: %.2f',rmse);
text(datenum(1993,5,1),3.2,rmse,'FontSize',15)
bias = nanmean(W(1).hs_ts(b2) - D.hs(b2));
bias = sprintf('Bias: %.2f',bias);
text(datenum(1993,5,1),2.8,bias,'FontSize',15);
title('Wind speeds: 5 - 10 m/s')

% 10 - 15 m/s bin
subplot(4,1,3)
o = plot(D.time_i(b3),D.hs(b3));
datetick()
hold on 
p = plot(W.time(b3),W.hs_ts(b3),'o','color',[.7 .7 .7]);
% lgd = legend([o,p],'Sentry Shoal Obs','LUT Predictions','Location','NorthEast');
% lgd.FontSize = 12;
ylabel('Wave Height [m]')
set(gca,'FontSize',10)
grid on
rmse = sqrt(nanmean((W.hs_ts(b3) - D.hs(b3)).^2));
rmse = sprintf('RMSE: %.2f',rmse);
text(datenum(1993,5,1),5.0,rmse,'FontSize',15)
bias = nanmean(W(1).hs_ts(b3) - D.hs(b3));
bias = sprintf('Bias: %.2f',bias);
text(datenum(1993,5,1),4.4,bias,'FontSize',15);
title('Wind speeds: 10 - 15 m/s')


% 10 - 15 m/s bin
subplot(4,1,4)
o = plot(D.time_i(b4),D.hs(b4));
datetick()
hold on 
p = plot(W.time(b4),W.hs_ts(b4),'o','color',[.7 .7 .7]);
% lgd = legend([o,p],'Sentry Shoal Obs','LUT Predictions','Location','NorthEast');
% lgd.FontSize = 12;
ylabel('Wave Height [m]')
set(gca,'FontSize',10)
grid on
rmse = sqrt(nanmean((W.hs_ts(b4) - D.hs(b4)).^2));
rmse = sprintf('RMSE: %.2f',rmse);
text(datenum(1993,5,1),4.8,rmse,'FontSize',15)
bias = nanmean(W(1).hs_ts(b4) - D.hs(b4));
bias = sprintf('Bias: %.2f',bias);
text(datenum(1993,5,1),4.1,bias,'FontSize',15);
title('Wind speeds: > 15 m/s')

%% QQ Scatter Analysis 
clf
cdf1 = sort(D.hs,'ascend');
cdf1y = linspace(0,1,length(cdf1));

cdf2 = sort(W.hs_ts,'ascend');
cdf2y = linspace(0,1,length(cdf2));

%Interp cdfs
cdf2i = interp1(cdf2y,cdf2,cdf1y);

x_axis = linspace(0,1,length(cdf1));
scatter(D.hs,W.hs_ts)
hold on
plot(cdf1,cdf2i,'r+')
grid on
hold on
line([0 6],[0 6],'Color','k','LineStyle','--')
xlabel('Observed Wave Heights [m]')
ylabel('Modeled Wave Heights [m]')
set(gca,'FontSize',14)



%% Target Plots 

% -------------- calculate target stats
hs = cell(length(W),1);
hs = cell(1,1);
hs{1} = cell(length(W),1);
for ii = 1:length(W)
    hs{1}{ii} = W(ii).hs_ts;
end
obs = cell(1,1);
obs{1} = H.Hs;

% dmwl=cellfun(@mean,wli); % Mean of the OBS
om=mean(H.Hs); % Mean of OBS
mhs=cellfun(@(x)(cell2mat(x')),hs,'un',0); % Converts the cell array for the station water levels to a 80000x6 matrix.  Turns the cell into a matrix
mmhs=cellfun(@(x)(nanmean(x,1)'),mhs,'un',0); % Take the mean of each row of the new matrix - Mean of the PREDICTED 

bias=cellfun(@(x)(x-om),mmhs,'un',0)'; % Calculates Bias - PREDICTED - OBSERVED 

temp_dhsu=cellfun(@(x,y)(x-y),obs,num2cell(om),'un',0); %remove mean from obs
%- This is necessary to make it work 
temp_dhsu = temp_dhsu{1,1}; 
dhsu = cell(length(W),1);
for ii = 1:length(W)
    dhsu{ii} = temp_dhsu;
end
temp_obs = cell(length(W),1);
for ii = 1:length(W)
    temp_obs{ii} = H.Hs;
end
obs = temp_obs;


crmsd=cell(length(bias),1);
rmsd=cell(length(bias),1);
for i=1:length(bias)
    mhsu=num2cell(bsxfun(@minus,mhs{i},mmhs{i}'),1)'; %remove mean from mod
    crmsd{i}=cellfun(@(x,y)(sqrt((sum((x-y).^2))/numel(x))),mhsu,dhsu).*... %craziness
        cellfun(@(x,y)(sign(std(x)-std(y))),num2cell(mhs{i},1)',obs); 
    rmsd{i}=cellfun(@(x,y)(sqrt((sum((x-y).^2))/numel(x))),...
        num2cell(mhs{i},1)',obs)';
end

% ------------  plot the target diagram
crmsdm=cell2mat(crmsd');
biasm=cell2mat(bias');
maxr=max([crmsdm(:);biasm(:)]);
maxr2=maxr+0.05*maxr;

xi=linspace(0,2*pi,100);
r=0:0.1:0.5;
[x,y]=cellfun(@(x)(pol2cart(xi,x)),num2cell(r),'un',0);

f=figure;
set(f,'renderer','zbuffer','units','inches',...
    'position',[1 1 8 7],...
    'paperpositionmode','au',...
    'color','w','inverthard','off')

hold on 
cellfun(@(x,y)(plot(x,y,'k--')),x,y)

plot(crmsdm',biasm','color',[0.6 0.6 0.6])

% -------------- Add Colored Dots to Target Plot

bc = parula(length(W));
for i=1:length(W)
    dh(i)=plot(999,999,'o',...
        'color',bc(i,:),'markerfacecolor',bc(i,:),...
        'linestyle','none');
end
% Andrew plots these way away so that he can call them in the legend I
% think 

cols=num2cell(parula(length(W)),2);
for i=1:length(W)
   lh(i)=plot(crmsdm(i),biasm(i),'o',...
        'color',bc(i,:),'markerfacecolor',bc(i,:),'linestyle','none');
end


% ---------------- Add Legend

lh2=lh(1,:);
lht=[dh(:);lh2(:)];
legstr={'H1';'H2';'H3';'H4';'H5';'H6';'H7';'H8';'H9';};

leg=legend(lht,legstr);
set(leg,'location','northeastoutside',...
    'fontang','it','fontweight','b',...
    'box','off','interpreter','none','FontSize',14)

set(gca,'da',[1 1 1],...
    'xaxislocation','origin',...
    'yaxislocation','origin',...
    'xlim',[-0.6 0.6],...
    'ylim',[-0.6 0.6],...
    'xtick',r,'ytick',r)

xl=xlim;
text(-.4,0,'\bf\itRMSD''(m)',...
    'horizontalalign','right',...
    'verticalalign','bot')
text(0,xl(2),'\bf\itBias (m)',...
    'horizontalalign','right',...
    'rotation',90,...
    'verticalalign','top')

