function [A,rx,ry] = Whidbey_RectiMerge(file1,file2)
%  I need a function that rectifies 2 cam views and blends them, given a
%  time and image type from the Argus products.
%
% inputs:   time or filename?... filename will be good. 
%           vertical level to rectify onto (e.g. tide)
%           save image or not? 
%
% outputs: rectified mosiac
%
% For the moment, let's keep this function specific to SunsetArgus.... 
% which will save all the inputs, etc. 
% 
% Inputs needed for rectification:
% gcps
% geom solution  % These "betas" are theoretically determined on the
% reference images and (should) not change.... Maybe just load the betas (they live in meta.betas)!
% srhREctImage wants time, images x,y,z, betas(1), globs (for lcp and knowns)
%inputs.rectxy = [50 0.5 1150 -1000 0.5 100];
%inputs.rectz = %tide or 0 or user defined

% file1 = 'E:\Puge\groundcontrol\1526763601.c1.timex.jpg'
% file2 = 'E:\Puge\groundcontrol\1526763601.c2.timex.jpg'
% file1 = 'E:\Puge\geomDay2\1526767203.c1.bright.jpg'
% file2 = 'E:\Puge\geomDay2\1526767203.c2.bright.jpg'

% 
% file1 = 'I:\Argus\Whidbey_2018\groundControl\1526763601.c1.timex.jpg'
% file2 = 'I:\Argus\Whidbey_2018\groundControl\1526763601.c2.timex.jpg'

  I1 = imread(file1); 
  I2 = imread(file2); 
% Meta data from previously established geometry of Whidbey c1 and c2:
  Mc1 = load(['I:\Argus\Whidbey_2018\groundControl\1526763601.c1.timex.Whidbey_c1.meta.mat']);
  Mc2 = load(['I:\Argus\Whidbey_2018\groundControl\1526763601.c2.timex.Whidbey_c2.meta.mat']);  
  
%%  Rectify the images:
  rectxy = [0 0.1 600 -700 0.1 0]; % [x0 dx x1 y0 dy y1] of rectified domain  "FUll View"
%   rectxy = [0 0.1 300 -450 0.1 5]; % [x0 dx x1 y0 dy y1] of rectified domain  "mid View"
%   rectxy = [0 0.05 70 -150 0.05 -50]; % [x0 dx x1 y0 dy y1] of rectified domain  "tight beach View"
  
%   if nargin < 3;     % vertical rect level
    rectz = 0; % vertical level onto which to rectify  
%   end
  

 % cam1:
  beta1 = Mc1.meta.betas;  
  lcp1 = Mc1.meta.globals.lcp;
  %known = Mc1.meta.globals.knownFlags
  [rectI1,rx,ry,rxyz] = srhRect(rectxy, rectz, I1, beta1, lcp1);

 % cam2:
  beta2 = Mc2.meta.betas;  
  lcp2 = Mc2.meta.globals.lcp;
  [rectI2,~,~] = srhRect(rectxy, rectz, I2, beta2, lcp2);
 

% figure; imshow(rectI1)
% figure; imshow(rectI2)

%% Merge the two retified camera views;
% 
%  Currently, this uses my super strange mask method, rather than the slick
%  Argus way... Eventually I'd like to figure outhow the argus function
%  works and implement it instead.
  [ni,mi,ci]=size(rectI1);
%    COMBO = nan(size(IPc1.finalImages.rectI)); 
%    COMBOy = IPc1.finalImages.y; COMBOx = IPc1.finalImages.x;
%    diff(COMBOy)
 % FOV boundaries!! (usually a fit object) 
 % X&Y points on left and right lines:
 
 % % % Try in array indices:
%   % "FULL VIEW":
% pts = [750 7001; 6001 6649 ; 247 6530; 4474 1; 327 6566; 4760 1; 1 1; 1 6483]';  % This is for rectxy = [0 0.1 600 -700 0.1 0];   
   pts = [750 7001; 6001 6649 ; 145 6686; 4474 1; 327 6566; 4760 1; 1 1; 1 6483]';  % This is for rectxy = [0 0.1 600 -700 0.1 0];
%   % "MID VIEW":
%    pts = [351 4535; 3001 4358; 143 4183; 2853 1;...
%           327 4068; 2991 129; 1 3983; 1 1]';  % This is for rectxy = [0 0.1 300 -450 0.1 5];
%   % "Tight Beach IEW":
%    pts = [530 2001; 1401 2001; 530 2000; 1401 671;...
%           742 2001; 1401 1015; 1 2001; 1 1]';  % This is for rectxy = [0 0.05 70 -150 0.05 -50];
      
      
      
 % pts = [57 2033; 2111 2201 ; 1 1864; 1567 1; 65 1827; 1727 1; 1 503; 38 1]';
 %Lines; % assume all y's in the domain, secribe Lines by x = (y-b)/m
 % But, don't forget that x&y are cordindates on the figure of the local
 % domain, and not local coordinate x&y's!! So, COMBOy is actually 'x' on
 % the figure (alongshore on the bottom axis)!! 
  y = 1:1:ni;
  x = 1:1:mi;
 %  y = COMBOx; 
 %  x = COMBOy; 
 for li = 1:4
   m(li) = (pts(2,2*li)-pts(2,2*li-1))/(pts(1,2*li)-pts(1,2*li-1)); 
   xA{li} = (y'-pts(2,2*li-1))./m(li) + pts(1,2*li-1);
   yA{li} = m(li).*(x' - pts(1,2*li-1)) + pts(2,2*li-1);
   b{li} = m(li) * (0 - pts(1,2*li)) + pts(2,2*li);   
 end
 %Image Mask the non-FOV space:
  I1 = rectI1; I1_dub = double(I1); 
  I2 = rectI2; I2_dub = double(I2); 
  [ni,mi,ci]=size(I1);
    mask_C1 = true(ni,mi); size(mask_C1);
    for ix = 1:length(x)  % Cam1 Left Boundary  
        iy = find(y >= yA{1}(ix));
      mask_C1(iy,ix) = false;
    end  
%     for iy = 1:length(y)  % Cam1 Left Boundary  
%         ix = find(x <= xA{1}(iy));
%       mask_C1(iy,ix) = false;
%     end  
    for iy = 1:length(y)  % Cam1 Right Boundary
      ix = find(x <= xA{2}(iy));
      mask_C1(iy,ix) = false;
    end  
%     im1 = imagesc(I1)
%     hold on; plot(xA{5},xA{6},'.r')
%     im1.AlphaData = mask_C1; 
    
    mask_C2 = true(ni,mi); size(mask_C2);
    for iy = 1:length(y)  % Cam2 Left Boundary  
      ix = find(x >= xA{3}(iy));
      mask_C2(iy,ix:end) = false;
    end  
    for iy = 1:length(y)  % Cam2 Right Boundary
      ix = find(x <= xA{4}(iy));
      mask_C2(iy,ix) = false;
    end      
%     im2.AlphaData = mask_C2; 
    %expand masks to images dimensions     
      for ci = 1:3
        maskC1(:,:,ci) = mask_C1(:,:,1);
        maskC2(:,:,ci) = mask_C2(:,:,1);
      end     
      %NaN images values where mask is true
      I1_dub(~maskC1) = NaN;
      I2_dub(~maskC2) = NaN;    
    %%% VECTORIZE GRIDS AND WORK OUT INTERSECTING AREA 
    %%% RE-ARRANGE EACH GRIDS ONTO ROW VECTORS
     vec_C1 = I1_dub(:)';
     vec_C2 = I2_dub(:)';
    %%%% LOCATE AREA ON THE EXTENDED GRIDS WHERE THE TWO IMAGES OVERLAY
     id_intersect = ~isnan(vec_C1) & ~isnan(vec_C2);
     id_C1 = ~isnan(vec_C1);
     id_C2 = ~isnan(vec_C2);
    % Take 2D only:
     id_intersect_2d = id_intersect(1:ni*mi);

    %%%% MESHGRID [X,Y] TO FIND OUT EACH POINT'S COORDINATES ON GRIDS
     [xx,yy] = meshgrid(x,y);
     vec_X = xx(:);
     vec_Y = yy(:);

   %%%% WEIGHT OVERLAYED AREA ACCORDING TO DISTANCES TO EACH IMAGE'S BOARDER
   % Lines   [slope, y-int]  (like the coefficients of polyval's p)
   % Convert line coefficients to the ax+by+c=0 form: 
    % y = Mx + B   --> Mx - y +B = 0
    % ax + by + c = 0      
    LINE_rb_C1 = [m(2),  b{2} ] ; % [slope, y-int]
    LINE_lb_C2 = [m(3), b{3} ] ;
    %  Dist = myPnt2LineDist(Pnt,Line,WhichAxis)  
    % Points
    xp = vec_X(id_intersect_2d);
    yp = vec_Y(id_intersect_2d);
    PNT = [xp yp];
    % Distances from each point within overlayed area and images' boarders
    % Check the distance between every point in the intersection of Cam1 and Cam2...
    dist_rb_C1 = myPnt2LineDist(PNT,LINE_rb_C1,'x'); 
    dist_lb_C2 = myPnt2LineDist(PNT,LINE_lb_C2,'x');
    dist_tot = dist_rb_C1+dist_lb_C2;
    % Attribute weights
    w_rb_C1 = dist_rb_C1./dist_tot;
    w_rb_C1_3d = repmat(w_rb_C1',1,3);
    w_lb_C2 = dist_lb_C2./dist_tot;
    w_lb_C2_3d = repmat(w_lb_C2',1,3);
    % Expand vectors with weights so it will end up equivalent to 3D vector 
    % (For plotting)
    vec_w_rb_C1 = nan(1,ni*mi);
    vec_w_rb_C1(id_intersect_2d) = w_rb_C1;
    vec_w_rb_C1 = repmat(vec_w_rb_C1,1,3);
    vec_w_lb_C2 = nan(1,ni*mi);
    vec_w_lb_C2(id_intersect_2d) = w_lb_C2;
    vec_w_lb_C2 = repmat(vec_w_lb_C2,1,3);
    % Calculate weighted averages over overlayed area
    vec_overlayC1 = vec_C1(id_intersect);
    vec_overlayC2 = vec_C2(id_intersect);
    vec_overlay_weightAveraged = (vec_overlayC1.*w_rb_C1_3d + vec_overlayC2.*w_lb_C2_3d); 
    %%% DEFINE AND RE-ARRANGE FINAL IMAGE
    % Final vectorized image
    vec = nanmean([vec_C1; vec_C2]);
    vec(id_intersect) = vec_overlay_weightAveraged;
    % To plot cameras masks:
    vec1 = vec_C1;
    vec1(~isnan(vec_C1)) = 1;
    vec1(isnan(vec_C1)) = 0;
    vec1(id_intersect) = w_rb_C1_3d;
    vec2 = vec_C2;
    vec2(~isnan(vec_C2)) = 1;
    vec2(isnan(vec_C2)) = 0;
    vec2(id_intersect) = w_lb_C2_3d;
    % %locate and replace NaNs for desirable gray scale.
     grayscale = 0; %0 is black background, 255 is white background
     inan = isnan(vec);
     vec_noNan = vec;
     vec_noNan(inan) = 1*grayscale;
    % %re-arrange vectors into images
    % A = uint8(reshape(vec,ny,nx,nc)); % FINAL MOSAIC
    A = uint8(reshape(vec_noNan,ni,mi,ci)); % FINAL MOSAIC
end
 %%
%  
% %  
% %  figure(13), set(13, 'Pos',[3200 200 560 420].*[1 1 2 2]); clf
% % % im2 = imagesc(A);
% figure; 
% % im2 = imagesc(ry,rx,rot90(A));
% im2 = imagesc(rx,ry,A);
% axis xy; axis image; grid off;
% view(-90,90)
% % im2 = imagesc(ry,rx,flipud(rot90(A)));
% ylabel('alongshore (m)'); xlabel('cross-shore (m)');
% %xlabel('alongshore (m)'); ylabel('cross-shore (m)');
% % title([Mc2.meta.globals.lcp.station ' merged timex: ' datestr(IPc2.finalImages.dn)])
% % axis xy; axis image; grid off;
% set(gca,'XDir','no');
% set(gca,'TickDir','out')
% set(gca,'XTick',[0:20:600],'Ytick',[-700:20:100])
% set(gca,'YTickLabels',abs(get(gca,'YTick')))
% % % Full view:
% % xlim([0 500])   
% % ylim([-600 00])
% % % Mid view:
% xlim([0 300])   
% ylim([-450 5])
% % % Tight Beach view:
% % xlim([0 70])   
% % ylim([-150 -50])
% 
% 
% % num2str(ylim)
% title('Whidbey 201805191400')
%  
%  save('tightBeachView.mat','rx','ry','A'); 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 