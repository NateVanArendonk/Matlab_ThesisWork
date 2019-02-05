


%%  Open the Stack, sort through it's instruments, plot the instruments.

%1540152000.c2.stack.ras
%1540152000.c1.stack.ras
stackout = 'I:\Argus\Whidbey_2018\data\timestacks\png\'

%% Here's how to explore the content of a single timestack file, e.g. c2:
  % StackFileName:
    filename = fullfile('I:\Argus\Whidbey_2018\data\timestacks','1540152000.c2.stack.ras')    
    % This spits out the params, time of each frame, frame number, and data
     % for each UV.... data(1,9) is 9th UV, 1st entry. Still, there is no
     % definitino of the UV positions to decode the stack.... I think we need
     % to use the insts list. 
      %[c1.p, c1.epoch, c1.MSC, c1.data] = loadStack( stack1 ); size(c1.data);   %  These stacks are especially large - too far offshore maybe?! 
     [c2.p, c2.epoch, c2.MSC, c2.data] = loadStack(filename); size(c2.data)
     hold on; 
     puvs = plot(c2.p.U,c2.p.V,'or')
     
     % Hunt for the lines using diff thresholding:
     %figure; plot(diff(c2.p.U))
       unique(find(diff(c2.p.U)>1))
       dinds = find(diff(c2.p.U)>1)
       %c2.p.U(1:dinds(1))
     for ii = 1:length(dinds)+1
      if ii == 1
       plot(c2.p.U(1:dinds(ii)),c2.p.V(1:dinds(ii)),'.k')
      elseif ii == length(dinds)+1
        plot(c2.p.U(dinds(ii-1)+1:end),c2.p.V(dinds(ii-1)+1:end),'.k')
      else
        plot(c2.p.U(dinds(ii-1)+1:dinds(ii)),c2.p.V(dinds(ii-1)+1:dinds(ii)),'.k')
      end
      pause
     end   % Looks ok... 
     % Define the lines
     
     
     
% Actually, the lines data is saved in the variable 'r' used to define the
% pixel lists... Let's use those to avoid repeated work.
   load('I:\Argus\Whidbey_2018\fieldTools\whidbeyPIX.c2.r.mat')  % r
   figure; plot(r.cams.XYZ(:,1),r.cams.XYZ(:,2),'or')
     %figure; plot(diff(r.cams.XYZ(:,1)))
     dinds = find(diff(r.cams.XYZ(:,1))<-5)
     hold on
     for ii = 1:length(dinds)+1
      if ii == 1
       plot(r.cams.XYZ(1:dinds(ii),1),r.cams.XYZ(1:dinds(ii),2),'.k')
      elseif ii == length(dinds)+1
        plot(r.cams.XYZ(dinds(ii-1)+1:end,1),r.cams.XYZ(dinds(ii-1)+1:end,2),'.k')
      else
        plot(r.cams.XYZ(dinds(ii-1)+1:dinds(ii),1),r.cams.XYZ(dinds(ii-1)+1:dinds(ii),2),'.k')
      end
      pause
     end   % Looks ok... 
     
     % Does it hold for UVs?
     figure; plot(r.cams.U,r.cams.V,'or')
     %figure; plot(diff(r.cams.XYZ(:,1)))
     dinds = find(diff(r.cams.XYZ(:,1))<-5);
     hold on
     for ii = 1:length(dinds)+1
      if ii == 1
       plot(r.cams.U(1:dinds(ii)),r.cams.V(1:dinds(ii)),'.k')
      elseif ii == length(dinds)+1
        plot(r.cams.U(dinds(ii-1)+1:end),r.cams.V(dinds(ii-1)+1:end),'.k')
      else
        plot(r.cams.U(dinds(ii-1)+1:dinds(ii)),r.cams.V(dinds(ii-1)+1:dinds(ii)),'.k')
      end
      pause
     end   % Looks ok... 
     clear inds
     for ii = 1:length(dinds)+1
         if ii == 1
           inds(ii,:) = [1 dinds(ii)];
         elseif ii == length(dinds)+1
           inds(ii,:) = [dinds(ii-1)+1 length(r.cams.U)];             
         else
           inds(ii,:) = [dinds(ii-1)+1 dinds(ii)];
         end
     end
     % inds   [startIndice endIndice] of each line
     %      % We should really pad the c2 variable with info from r. Then c1 and
%      % c2 can be easily differentiated.  I'm not sure if r.cams.U and V are
%      % different that c2.p.U and V? They might be?
%      
%      size(c2.p.U)
%      size(r.cams.U)
%      % Hmm.. 1 less point in the r veector... I thin kI recall there being
%      % an extra header line or something in the p.      
%      size(x)
%      size(t)
%      size(c2.data(:,inds(ii,1):inds(ii,2)))
%      %  Well, it doesn't really matter if you cruise through the positions
%      %  to find which are the lines.
        c2.r.inds = inds; 
        c2.r.X = r.cams.XYZ(:,1);
        c2.r.Y = r.cams.XYZ(:,2); 
        c2.r.U = r.cams.U;
        c2.r.V = r.cams.V;
     
     % E.g. for line 7, ii = 7
     for ii = 1:length(c2.r.inds)
     %figure; imshow(c2.data(:,inds(ii,1):inds(ii,2)))  % In uv coords.
     size(c2.data(:,c2.r.inds(ii,1):c2.r.inds(ii,2)))
     t= c2.epoch-c2.epoch(1); % relative time in seconds
     datetime0 = datetime(datetime(datestr(epoch2Matlab(c2.epoch(1))),'TimeZone','UTC'),'TimeZone','America/Los_Angeles')
     ltime0 = datenum(datetime0); datestr(ltime0)  % Local time
     x = c2.r.X(c2.r.inds(ii,1):c2.r.inds(ii,2));
     fig = figure; pcolor(x,t, c2.data(:,c2.r.inds(ii,1):c2.r.inds(ii,2)))  % In uv coords.
     set(fig,'InvertHardCopy','off')
     set(fig,'Color','k')
     set(gca,'XColor','w','yColor','w','Layer','top')
     
 
     shading interp; colormap(gray)
     axis ij 
     caxis([40 240])
     xlabel('x (m)')
     ylabel(['t (s) after ' datestr(ltime0,'dd-mmm-yyyy HH:MM') 'PDT'])
     title(['x-shore transect ' num2str(ii)],'color','w')
     set(gcf,'Pos',[3300 100 600 1400])
     drawnow 
     print('-dpng','-r150',[stackout 'whidbey.tran' num2str(ii,'%02.0f') '.' num2str(c2.epoch(1),'%10.0f') '.c2.png'  ])
     pause(.1)
     ylim([0 180])
     print('-dpng','-r150',[stackout 'whidbey.tran' num2str(ii,'%02.0f') '.' num2str(c2.epoch(1),'%10.0f') '.c2.ZOOM.png'  ])
     close gcf
     end
     

     
     
     %%  Now let's look at a c1 file, and think about how to join the lines that span both cameras FOVs. 
     
     % c1 only. 
       filename = fullfile('I:\Argus\Whidbey_2018\data\timestacks','1540152000.c1.stack.ras')    
       [c1.p, c1.epoch, c1.MSC, c1.data] = loadStack(filename); size(c1.data)
       load('I:\Argus\Whidbey_2018\fieldTools\whidbeyPIX.c1.r.mat')  % r
        dinds = find(diff(r.cams.XYZ(:,1))<-5);
        clear inds
        for ii = 1:length(dinds)+1
         if ii == 1
           inds(ii,:) = [1 dinds(ii)];
         elseif ii == length(dinds)+1
           inds(ii,:) = [dinds(ii-1)+1 length(r.cams.U)];             
         else
           inds(ii,:) = [dinds(ii-1)+1 dinds(ii)];
         end
        end
        c1.r.inds = inds; 
        c1.r.X = r.cams.XYZ(:,1);
        c1.r.Y = r.cams.XYZ(:,2); 
        c1.r.U = r.cams.U;
        c1.r.V = r.cams.V;
     
       for ii = 1:length(c1.r.inds)        
         size(c1.data(:,c1.r.inds(ii,1):c1.r.inds(ii,2)))
         t= c1.epoch-c1.epoch(1); % relative time in seconds
         datetime0 = datetime(datetime(datestr(epoch2Matlab(c1.epoch(1))),'TimeZone','UTC'),'TimeZone','America/Los_Angeles')
         ltime0 = datenum(datetime0); datestr(ltime0)  % Local time
         x = c1.r.X(c1.r.inds(ii,1):c1.r.inds(ii,2));
         figure; pcolor(x,t, c1.data(:,c1.r.inds(ii,1):c1.r.inds(ii,2)))  % In uv coords.
         shading interp; colormap(gray)
         axis ij 
         caxis([40 240])
         xlabel('x (m)')
         ylabel(['t (s) after ' datestr(ltime0,'dd-mmm-yyyy HH:MM') 'PDT'])
         title(['x-shore transect ' num2str(ii)])
         set(gcf,'Pos',[3300 100 600 1400])
         drawnow 
         print('-dpng','-r150',[stackout 'whidbey.tran' num2str(ii,'%02.0f') '.' num2str(c1.epoch(1),'%10.0f') '.c1.png'  ])
         pause(.1)
         ylim([0 180])
         print('-dpng','-r150',[stackout 'whidbey.tran' num2str(ii,'%02.0f') '.' num2str(c1.epoch(1),'%10.0f') '.c1.ZOOM.png'  ])
         close gcf
       end
     
     
     
       %% Now how about joining the partial lines?
       % line 5 and line 6 (c1) are also line 1 and line 2 (c2)
       
    % Transect 5   
    xc = [c2.r.X(c2.r.inds(1,1):c2.r.inds(1,2))  ; c1.r.X(c1.r.inds(5,1):c1.r.inds(5,2))];
    Ic = [c2.data(:,c2.r.inds(1,1):c2.r.inds(1,2)) c1.data(:,c1.r.inds(5,1):c1.r.inds(5,2))];
    
     figure; pcolor(xc,t, Ic)  % In uv coords.
         shading interp; colormap(gray)
         axis ij 
         caxis([40 240])
     %overlap region:
     line([ dup(max(c2.r.X(c2.r.inds(1,1):c2.r.inds(1,2))),2)], ylim, 'color','y')
     line([ dup(min(c1.r.X(c1.r.inds(5,1):c1.r.inds(5,2))),2)], ylim, 'color','r')
     
    % Transect 5
    xc = [c2.r.X(c2.r.inds(2,1):c2.r.inds(2,2))  ; c1.r.X(c1.r.inds(6,1):c1.r.inds(6,2))];
    Ic = [c2.data(:,c2.r.inds(2,1):c2.r.inds(2,2)) c1.data(:,c1.r.inds(6,1):c1.r.inds(6,2))];
      figure; pcolor(xc,t, Ic)  % In uv coords.
         shading interp; colormap(gray)
         axis ij 
         caxis([40 240])
     %overlap region:
     line([ dup(max(c2.r.X(c2.r.inds(2,1):c2.r.inds(2,2))),2)], ylim, 'color','y')
     line([ dup(min(c1.r.X(c1.r.inds(6,1):c1.r.inds(6,2))),2)], ylim, 'color','r')
    


%%  Let's now plot a rectifed-merged image product from the same time with and run-up line for reference:

  file1 = fullfile('I:\Argus\Whidbey_2018\data\products','1540152000.c1.snap.jpg');
  file2 = fullfile('I:\Argus\Whidbey_2018\data\products','1540152000.c2.snap.jpg');
  rectz = 0;   
    estart = regexp(file1,'\d{10}')    
    producttype = file1(estart+14: regexp(file1,'.jpg')-1);
    epoch = str2num(file1(estart:estart+9))
    localdatetime = datetime(datetime(datestr(epoch2Matlab(epoch)),'timezone','UTC'),'timezone','America/Los_Angeles');
    localtime = datenum(localdatetime); datestr(localtime) 
  [A,rx,ry] = Whidbey_RectiMerge(file1,file2);        
  [fig] = plotWhidbeyRectMergeBLACK(A,rx,ry,rectz,epoch2Matlab(epoch)); 
  set(fig, 'Pos',[3200 120 1366 710].*[1 1 1 1]);
   xlim([-300 0])   
   ylim([0 200])
   set(gca,'XTick',[-400:50:50],'Ytick',[0 :50:200])
   set(gca,'XTickLabels',abs(get(gca,'XTick')));
   hold on; 

   
   k=0;
   if exist('pline'); delete(pline); end
   for ci = 1:2
     for ii = 1:length(eval(['c' num2str(ci) '.r.inds']))                    
          X = eval(['c' num2str(ci) '.r.X(c' num2str(ci) '.r.inds(ii,1): c' num2str(ci) '.r.inds(ii,2))']);
          Y = eval(['c' num2str(ci) '.r.Y(c' num2str(ci) '.r.inds(ii,1): c' num2str(ci) '.r.inds(ii,2))']);
          k = k+1;
          pline(k) = plot(Y,X,'.y');          
     end
   end
       
   print('-dpng','-r150',[stackout 'whidbey.tranMap.' num2str(epoch,'%10.0f') '.cx.' producttype '.png'  ])
        
      
   
   
   
