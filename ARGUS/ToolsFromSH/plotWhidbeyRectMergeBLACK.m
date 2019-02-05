function [fig] = plotWhidbeyRectMergeBLACK(A,rx,ry,rz,dnGMT)


%  SBm.TIME(SBinds), SBm.ECg(SBinds), SBm.Dp(SBinds), SBm.Sxx(SBinds), SBm.Sxy(SBinds), wlMont.time(MBinds), wlMont.WL_VALUE(MBinds)
% The figure: 
fig = figure; 
% set(fig, 'Pos',[-3200 120 830 750].*[1 1 1 1]); 
set(fig, 'Pos',[3200 120 1366 710].*[1 1 1 1]); 
set(fig,'InvertHardCopy','off')
% TIMEZONE STUFF! 
 % Define their timezone and convert to local
dtime = datetime(datetime(datestr(dnGMT),'TimeZone','UTC'),'TimeZone','America/Los_Angeles');  mdtime = datenum(dtime);% datestr(mdtime)
% SBmTIME =  datetime(datetime(datestr(metdata.SBmTIME),'TimeZone','UTC'),'TimeZone','America/Anchorage');
% wlMonTime = datetime(datetime(datestr(metdata.wlMonTime),'TimeZone','UTC'),'TimeZone','America/Anchorage');



% The Argus imagery: 
 ax(1) = axes; 
 set(ax(1),'Pos',[0.08   0.09    0.88    0.87])
 axes(ax(1))
  im2 = imagesc(ry,rx,flipud(rot90(A)));
   xlabel('alongshore (m)'); ylabel('cross-shore (m)');          
   axis xy; axis image; grid off; set(gca,'XDir','rev'); set(gca,'TickDir','out')
%    xlim([-200 10])      
%    ylim([10 200])
%    xlim([-200 0])   
%    ylim([10 110])
% set(ax(1),'XTick',[-200:20:60],'Ytick',[0 :20:200])
xlim([-500 0])   
ylim([0 300])   
set(ax(1),'XTick',[-500:50:0],'Ytick',[0 :50:300])
    set(ax(1),'XTickLabels',abs(get(gca,'XTick')));
    set(fig,'Color','k')
    set(ax(1),'XColor','w','yColor','w')
    tits = title(['Whidbey Island: ' datestr(dtime,'dd-mmm-yyyy HH:MM') ' PDT'],'color','w');       
    
%    % ADd a compass for local coordinates:
%    ax(2) = axes; 
%    set(ax(2),'Pos',[.82 .489 .1 .1],'color','none','vis','off');
%    axes(ax(2))   ;
%     DIRS = mod([0 90 180 225 270]+120,360);
%     [u,v]=calc_uv(ones(size(DIRS)),DIRS); [ul,vl]=calc_uv(ones(size(DIRS))+.1,DIRS);
%      qv= quiver(zeros(size(u)),zeros(size(u)),u,v); axis image;
%      qv.Color = 'k'; qv.LineWidth = 2; qv.MaxHeadSize = 1; qv.AutoScaleFactor = .9; set(ax(2),'Vis','off');
%      XL = [-1.2 1.2]; xlim(XL); ylim(XL);
%      qlabs = {'N','E','S','SW','W'};
%      for ii = 1:length(qlabs)
%        qtxt(ii)=text(ul(ii),vl(ii),qlabs{ii},'FontSize',14,'horizontal','cent','fontweight','bold');
%      end
    
    
    
    
 
% % Wave and tide subplot axes:
%   % Wave ECg and Dp
%   ha(3) = axes('Position',[.065 .7 .225 .25],'color','k','Box','on','YColor','w','XColor','w'); 
%   ha(6) = axes('Position',[.065 .7 .225 .25],'color','none','Box','off','YColor','w','XColor','w','XAxisLocation','top','XTickLabel',[],'YAxisLocation','right'); 
%   % Radiatino Stress
%   ha(1) = axes('Position',[.38 .7 .225 .25],'color','k','Box','on','YColor','w','XColor','w'); 
%   ha(4) = axes('Position',[.38 .7 .225 .25],'color','none','XColor','w','XAxisLocation','top','XTickLabel',[],'YAxisLocation','right','YColor',rgb('gold'),'box','off');
%   % wlev / tide
%   ha(2) = axes('Position',[.70 .7 .225 .25],'color','k','Box','on','YColor',rgb('orangered'),'XColor','w'); 
%   ha(5) = axes('Position',[.70 .7 .225 .25],'color','none','XColor','w','XAxisLocation','top','XTickLabel',[],'YAxisLocation','right','YColor',rgb('DodgerBlue'),'box','off');
  %%

  
% Plot the subplots:
   % CDIP: Hs or ECg
%     axes(ha(1)); 
%      cla; hold on;            
%       % Hs:
% %       plot(SBm.TIME(SBinds),SBm.Hs(SBinds),'-w','MarkerSize',16,'linewidth',4)
% %       xlim([dnGMT-timeBuff dnGMT+timeBuff])
% %       ylim(ha(1),HsYL); set(ha(1),'YTIck',[0:1:HsYL(2)])    
% %       ln = line([dnGMT dnGMT],ylim,'color',rgb('firebrick'),'linewidth',5);
% %       uistack(ln,'bottom');
% %       ylabel(ha(1),'H_s (m)')
%       % ECg:
%       plot(SBmTIME,metdata.SBmECg./1E3,'-w','MarkerSize',16,'linewidth',4)
%       xlim([dtime-days(timeBuff) dtime+days(timeBuff)])
%       ylim auto%ylim(ha(1),metdata.ECgYL./1E3); set(ha(1),'YTIck',[0:round(metdata.ECgYL(2)/1E3/4):metdata.ECgYL(2)/1E3])    
%       ln1 = line([dtime dtime],ylim,'color',rgb('firebrick'),'linewidth',5);
%       uistack(ln1,'bottom');
%       ylabel(ha(1),'Ec_g (kJ m^{-1} s^{-1})')
%       title(ha(1),'Local Waves (CDIP)','color','w')
%    % CDIP: Dp
%     axes(ha(4)); cla; hold on; 
%       plot(SBmTIME,metdata.SBmDp,'linestyle','none','marker','.','color',rgb('gold'),'MarkerSize',16)           
%       xlim([dtime-days(timeBuff) dtime+days(timeBuff)])
%       ylim([215 280]); set(ha(4),'YTIck',[90:45:270],'YTickLabel',{'E','SE','S','SW','W'})
%       line(xlim,[241 241],'linewidth',1,'color',rgb('orange'))
%       ylDp = ylabel(ha(4),'D_p (dir from)');
%       ylDp.Units = 'normalized'; ylDp.Position = [1.0315 .5 0];
%       set(ha(4),'XTickLabel',[])
%       
%    % CDIP: Sxx and Sxy
%     axes(ha(2)); 
%      cla; hold on;            
%       % Sxx:
%       plot(SBmTIME,metdata.SBmSxx,'linestyle','-','color',rgb('orangered'),'linewidth',4)
%       xlim([dtime-days(timeBuff) dtime+days(timeBuff)])
%       ylim auto%ylim(ha(2),metdata.SxxYL); set(ha(2),'YTIck',[0:metdata.SxxYL(2)/3:metdata.SxxYL(2)])    
%       ln2 = line([dtime dtime],ylim,'color',rgb('firebrick'),'linewidth',5);
%       uistack(ln2,'bottom');
%       ylabel(ha(2),'S_{xx} (N/m)')
%       title(ha(2),'Wave Radiation Stress','color','w')
%       crosstext = text(.02,.96,'shoreward','units','normalized','color',rgb('orangered'),'fontsize',14,'fontweight','bold','horizontalalignment','left');
%    % CDIP: Sxy
%     axes(ha(5)); cla; hold on; 
%        plot(SBmTIME,metdata.SBmSxy,'linestyle','-','color',rgb('dodgerblue'),'linewidth',4)           
%       xlim([dtime-days(timeBuff) dtime+days(timeBuff)])
%       ylim auto; %ylim(ha(5),[-.025 .025]); set(ha(5),'YTIck',[-.025:0.01:.025])   
%      line(xlim,[0 0],'linewidth',1,'color',rgb('dodgerblue'))
%       ylabel(ha(5),'S_{xy} (N/m)')
%       set(ha(5),'XTickLabel',[])
%       lrtext = text(.98,.5,{'right';' ';'left'},'units','normalized','color',rgb('dodgerblue'),'fontsize',14,'fontweight','bold','horizontalalignment','right');
% 
%   % tide / wlev
%       wlevYL = [-.5 2.2];
%       MBdatumSTR = {'MHHW','MHW','MSL','MLW','MLLW','NAVD'};
%       MBdatumST = [2.657 2.443 1.893 1.364 1.031 .988];
%       MBdatumNAV = MBdatumST-MBdatumST(end);
  
  
%     axes(ha(3)); 
%      cla; hold on;            
%       % wlev:
%       plot(wlMonTime,metdata.wlMonWLEV,'linestyle','-','color',rgb('silver'),'linewidth',4)
%       xlim([dtime-days(timeBuff) dtime+days(timeBuff)])
%       ylim(ha(3),wlevYL); set(ha(3),'YTIck',[-.5:.5:2.5])                       
%       ylabel(ha(3),'wlev (m)')
%       title(ha(3),'Tide','color','w')
%       ln3 = line([dtime dtime],ylim,'color',rgb('firebrick'),'linewidth',5);    
%       for ii = 1: length(MBdatumNAV)-1  
%         lnNAVD(ii) = line(xlim,[MBdatumNAV(ii) MBdatumNAV(ii)],ylim,'color',rgb('silver'),'linewidth',1);
%         uistack(lnNAVD(ii),'bottom');
%       end
%       uistack(ln3,'bottom');      
%     axes(ha(6)); 
%      cla; hold on; 
%      xlim(ha(6),datenum(get(ha(3),'Xlim')))
%      ylim(ha(6),get(ha(3),'Ylim'))
%      set(ha(6),'YTick',fliplr(MBdatumNAV(1:end-1)),'YTickLabel',fliplr(MBdatumSTR(1:end-1)))
%    
%     %Dateticks:
%    
%       linkaxes(ha(1:5),'x')
      
      
      
      
      %%
