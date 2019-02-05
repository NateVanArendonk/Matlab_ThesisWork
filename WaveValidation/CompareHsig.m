% Load in wave hindcast output
clearvars


% Load in Hindcast results 
w = dir('Hindcast_Output/*.mat');
for ii = 1:length(w)
    W(ii) = load(['Hindcast_Output/' w(ii).name]);
end

% Load in Hein Bank Data
h = dir('E:\Abbas\Model_Met_Forcings\HeinBank\*.mat');
inds2get = 2:5;
H = struct;
H.time = [];
H.Hs = [];

% Load 2005 - 2008 waves 
for ii = 1:length(inds2get) % 2005 - 2008
    temp = load(['E:\Abbas\Model_Met_Forcings\HeinBank\' h(inds2get(ii)).name]);
    H.time = vertcat(H.time,temp.time);
    H.Hs = vertcat(H.Hs,temp.Hs);
end
clear h ii w inds2get




%%
%Get obs for same time period as model

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

%% Plot Timeseries 
% ---------------------------- First Group of 4 ---------------------------
clf
subplot(4,1,1)

o = plot(H.time,H.Hs);
datetick()
hold on 
p = plot(W(1).time,W(1).hs_ts,'o','color',[.7 .7 .7]);
lgd = legend([o,p],'Hein Bank Obs','LUT Predictions','Location','NorthEast');
lgd.FontSize = 12;
ylabel('Wave Height [m]')
set(gca,'FontSize',10)
grid on
text(datenum(2006,5,1),3.6,'NNRP Point: H1')
rmse = sqrt(nanmean((W(1).hs_ts - H.Hs).^2));
rmse = sprintf('RMSE: %.2f',rmse);
text(datenum(2006,5,1),3.2,rmse)
bias = nanmean(W(1).hs_ts - H.Hs);
bias = sprintf('Bias: %.2f',bias);
text(datenum(2006,5,1),2.8,bias);
xticks([datenum(2005,1,1) datenum(2006,1,1) datenum(2007,1,1) datenum(2008,1,1)]);
xticklabels({'2005','2006','2007','2008'})

subplot(4,1,2)
novInds = find(month(H.time) == 11);
plot(H.time,H.Hs)
datetick()
hold on
plot(W(2).time,W(2).hs_ts,'o','color',[.7 .7 .7])
ylabel('Wave Height [m]')
set(gca,'FontSize',10)
grid on
text(datenum(2006,5,1),3.6,'NNRP Point: H2')
rmse = sqrt(nanmean((W(2).hs_ts - H.Hs).^2));
rmse = sprintf('RMSE: %.2f',rmse);
text(datenum(2006,5,1),3.2,rmse)
bias = nanmean(W(2).hs_ts - H.Hs);
bias = sprintf('Bias: %.2f',bias);
text(datenum(2006,5,1),2.8,bias);
xticks([datenum(2005,1,1) datenum(2006,1,1) datenum(2007,1,1) datenum(2008,1,1)]);
xticklabels({'2005','2006','2007','2008'})

subplot(4,1,3)
novInds = find(month(H.time) == 11);
plot(H.time,H.Hs)
datetick()
hold on
plot(W(3).time,W(3).hs_ts,'o','color',[.7 .7 .7])
ylabel('Wave Height [m]')
set(gca,'FontSize',10)
grid on
text(datenum(2006,5,1),3.6,'NNRP Point: H3')
rmse = sqrt(nanmean((W(3).hs_ts - H.Hs).^2));
rmse = sprintf('RMSE: %.2f',rmse);
text(datenum(2006,5,1),3.2,rmse)
bias = nanmean(W(3).hs_ts - H.Hs);
bias = sprintf('Bias: %.2f',bias);
text(datenum(2006,5,1),2.8,bias);
xticks([datenum(2005,1,1) datenum(2006,1,1) datenum(2007,1,1) datenum(2008,1,1)]);
xticklabels({'2005','2006','2007','2008'})

subplot(4,1,4)
novInds = find(month(H.time) == 11);
plot(H.time,H.Hs)
datetick()
hold on
plot(W(4).time,W(4).hs_ts,'o','color',[.7 .7 .7])
ylabel('Wave Height [m]')
xlabel('Time')
set(gca,'FontSize',10)
grid on 
text(datenum(2006,5,1),3.6,'NNRP Point: H4')
rmse = sqrt(nanmean((W(4).hs_ts - H.Hs).^2));
rmse = sprintf('RMSE: %.2f',rmse);
text(datenum(2006,5,1),3.2,rmse)
bias = nanmean(W(4).hs_ts - H.Hs);
bias = sprintf('Bias: %.2f',bias);
text(datenum(2006,5,1),2.8,bias);
xticks([datenum(2005,1,1) datenum(2006,1,1) datenum(2007,1,1) datenum(2008,1,1)]);
xticklabels({'2005','2006','2007','2008'})

printFig(gcf,'HeinVsLUT_1_4',[8.5 11],'png',300)
%% ------------------------- Next Group of 5 -------------------------------

clf
subplot(5,1,1)
novInds = find(month(H.time) == 11);
o = plot(H.time,H.Hs);
datetick()
hold on 
p = plot(W(5).time,W(5).hs_ts,'o','color',[.7 .7 .7]);
lgd = legend([o,p],'Hein Bank Obs','LUT Predictions','Location','NorthWest');
lgd.FontSize = 12;
ylabel('Wave Height [m]')
set(gca,'FontSize',10)
grid on
text(datenum(2008,2,1),5.6,'NNRP Point: H5')
rmse = sqrt(nanmean((W(5).hs_ts - H.Hs).^2));
rmse = sprintf('RMSE: %.2f',rmse);
text(datenum(2008,2,1),5,rmse)
bias = nanmean(W(5).hs_ts - H.Hs);
bias = sprintf('Bias: %.2f',bias);
text(datenum(2008,2,1),4.4,bias);
xticks([datenum(2005,1,1) datenum(2006,1,1) datenum(2007,1,1) datenum(2008,1,1)]);
xticklabels({'2005','2006','2007','2008'})
ylim([0 6])

subplot(5,1,2)
novInds = find(month(H.time) == 11);
plot(H.time,H.Hs)
datetick()
hold on
plot(W(6).time,W(6).hs_ts,'o','color',[.7 .7 .7])
ylabel('Wave Height [m]')
set(gca,'FontSize',10)
grid on
text(datenum(2008,2,1),5.6,'NNRP Point: H6')
rmse = sqrt(nanmean((W(6).hs_ts - H.Hs).^2));
rmse = sprintf('RMSE: %.2f',rmse);
text(datenum(2008,2,1),5,rmse)
bias = nanmean(W(6).hs_ts - H.Hs);
bias = sprintf('Bias: %.2f',bias);
text(datenum(2008,2,1),4.4,bias);
xticks([datenum(2005,1,1) datenum(2006,1,1) datenum(2007,1,1) datenum(2008,1,1)]);
xticklabels({'2005','2006','2007','2008'})
ylim([0 6])


subplot(5,1,3)
novInds = find(month(H.time) == 11);
plot(H.time,H.Hs)
datetick()
hold on
plot(W(7).time,W(7).hs_ts,'o','color',[.7 .7 .7])
ylabel('Wave Height [m]')
set(gca,'FontSize',10)
grid on
text(datenum(2008,2,1),5.6,'NNRP Point: H7')
rmse = sqrt(nanmean((W(7).hs_ts - H.Hs).^2));
rmse = sprintf('RMSE: %.2f',rmse);
text(datenum(2008,2,1),5,rmse)
bias = nanmean(W(7).hs_ts - H.Hs);
bias = sprintf('Bias: %.2f',bias);
text(datenum(2008,2,1),4.4,bias);
xticks([datenum(2005,1,1) datenum(2006,1,1) datenum(2007,1,1) datenum(2008,1,1)]);
xticklabels({'2005','2006','2007','2008'})
ylim([0 6])

subplot(5,1,4)
novInds = find(month(H.time) == 11);
plot(H.time,H.Hs)
datetick()
hold on
plot(W(8).time,W(8).hs_ts,'o','color',[.7 .7 .7])
ylabel('Wave Height [m]')
set(gca,'FontSize',10)
grid on 
text(datenum(2008,2,1),5.6,'NNRP Point: H8')
rmse = sqrt(nanmean((W(8).hs_ts - H.Hs).^2));
rmse = sprintf('RMSE: %.2f',rmse);
text(datenum(2008,2,1),5,rmse)
bias = nanmean(W(8).hs_ts - H.Hs);
bias = sprintf('Bias: %.2f',bias);
text(datenum(2008,2,1),4.4,bias);
xticks([datenum(2005,1,1) datenum(2006,1,1) datenum(2007,1,1) datenum(2008,1,1)]);
xticklabels({'2005','2006','2007','2008'})
ylim([0 6])

subplot(5,1,5)
novInds = find(month(H.time) == 11);
plot(H.time,H.Hs)
datetick()
hold on
plot(W(9).time,W(9).hs_ts,'o','color',[.7 .7 .7])
ylabel('Wave Height [m]')
xlabel('Time')
set(gca,'FontSize',10)
grid on 
text(datenum(2008,2,1),5.6,'NNRP Point: H9')
rmse = sqrt(nanmean((W(9).hs_ts - H.Hs).^2));
rmse = sprintf('RMSE: %.2f',rmse);
text(datenum(2008,2,1),5,rmse)
bias = nanmean(W(9).hs_ts - H.Hs);
bias = sprintf('Bias: %.2f',bias);
text(datenum(2008,2,1),4.4,bias);
xticks([datenum(2005,1,1) datenum(2006,1,1) datenum(2007,1,1) datenum(2008,1,1)]);
xticklabels({'2005','2006','2007','2008'})
ylim([0 6])

printFig(gcf,'HeinVsLUT_5_9',[8.5 11],'png',300)
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

