
addpath 'I:\Argus\Sunset_2017\tools'
% Determine the existing imagery data to use:
imagepath = 'I:\Argus\Whidbey_2018\data\products\';
rectpath = 'I:\Argus\Whidbey_2018\data\rect\';

rectz = 0

% 1526851803
% 1527105603
% Specify a certain Epoch to rectimerge:
epochs = [1526851800, ...
1526851801, ...
1526851803, ...
1527105600, ...
1527105601, ...
1527105603];

for ii = 1:length(epochs)
estr = num2str(epochs(ii)); 
c1name = ls([imagepath estr '.c1.*.jpg'])
c2name = ls([imagepath estr '.c2.*.jpg'])
producttype = c1name(15:regexp(c1name,'.jpg')-1);

c1epochs = str2num([c1name(:,1:9) '0']);  % still in GMT!
c2epochs = str2num([c2name(:,1:9) '0']);  % still in GMT!
c1mdates = epoch2Matlab(c1epochs);  % still in GMT!
c1lmdates = c1mdates - 7/24; % now in local (PDT)

 z0 = 0;
   
  file1 = fullfile(imagepath, c1name)
  file2 = fullfile(imagepath, c2name)

  dnGMT = c1mdates; %datestr(dnGMT)
  dnLocal = datenum(datetime(datetime(datestr(dnGMT),'TimeZone','UTC'),'TimeZone','America/Los_Angeles'));
  disp(['Processing imagery for ' datestr(dnLocal) 'PDT...'])  ;
  [A,rx,ry] = Whidbey_RectiMerge(file1,file2);    
 
  [fig] = plotWhidbeyRectMergeBLACK(A,rx,ry,rectz,dnGMT);
  xlim([-420 0])   
  ylim([0 220]) 

  %
  if 1 % %Save PngS
      pngpath = [rectpath 'png\special\'];
      fn = [estr '.cx.whidbey.' producttype '.png'];
      eval(['print -dpng -r150 ' pngpath fn]);
  end
end