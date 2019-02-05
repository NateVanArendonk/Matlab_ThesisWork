addpath C:\Functions_Matlab\time
addpath C:\Functions_Matlab
addpath C:\Functions_Matlab\cbrewer
addpath C:\Functions_Matlab\circ_stats
addpath C:\Functions_Matlab\mapping\kml

load('MACA_Stats_Data_rcp85_MAxsGrabbed_1.mat')

clf
% Plotting Dimensions
p_left = .06;
p_right = .05;
p_top = .1;
p_bot = .06;
p_spacing = .01;
p_wid = (1-p_right-p_left-p_spacing)/3;
p_height = (1-p_top-p_bot-p_spacing)/2;

% Grab U AND V for Future and Historic 
temp_u_hist = nan_mean(u_max(:,:,years1),3);
temp_v_hist = nan_mean(v_max(:,:,years1),3);
wnddir_dif_matrixU = wrap2360(90-((180/pi)*atan2(temp_u_hist,temp_v_hist))+180); % Convert to Nautical 

temp_u_fut = nan_mean(u_max(:,:,years2),3);
temp_v_fut = nan_mean(v_max(:,:,years2),3);
wnddir_dif_matrixV = wrap2360(90-((180/pi)*atan2(temp_u_fut,temp_v_fut))+180); % Convert to Nautical 


% ------------- Plot U
row = 1;
col = 0;
a1 = axes('position',[p_left+col*(p_wid+p_spacing) p_bot+row*(p_height+p_spacing) p_wid p_height]);
pcolor(lon,lat,wnddir_dif_matrixU);
shading interp
hold on
plot(wa_lon, wa_lat,'k')
axis equal
xlim([min(lon(:)) -121.9])
ylim([min(lat(:)) 49.1])
caxis([200 340])
ylabel('Degrees Latitude')
set(gca,'XTickLabel',[])

% ------------- Plot V
col = 1;
a2 = axes('position',[p_left+col*(p_wid+p_spacing) p_bot+row*(p_height+p_spacing) p_wid p_height]);
colormap(parula);
pcolor(lon,lat,wnddir_dif_matrixV);
shading interp
hold on
plot(wa_lon, wa_lat,'k')
axis equal
xlim([min(lon(:)) -121.9])
ylim([min(lat(:)) 49.1])
caxis([200 340])
set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])

% ------------- Plot Difference
col = 2;
a3 = axes('position',[p_left+col*(p_wid+p_spacing) p_bot+row*(p_height+p_spacing) p_wid p_height]);
N_c = 100;
mycolors = flipud(cbrewer('div','RdBu',N_c));
cmap = colormap(a3,mycolors);
dif_mat = wnddir_dif_matrixV - wnddir_dif_matrixU;
pcolor(lon,lat,dif_mat);
shading interp
hold on
plot(wa_lon, wa_lat,'k')
axis equal
xlim([min(lon(:)) -121.9])
ylim([min(lat(:)) 49.1])
set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])
caxis([-10 10])

% % ------------- Add Colorbars
% cb1 = colorbar(a1,'Location','NorthOutSide','Position',[p_left .7 (2*p_wid+p_spacing) 0.025]);
% cb1.FontSize = 10;
% cb1.Label.String = 'Max Wind Direction [degrees]';
% cb1.Label.FontSize = 12;
% 
% cb2 = colorbar(a3,'Location','NorthOutSide','Position',[(p_left+2*p_wid+2*p_spacing) .7 p_wid 0.025]);
% cb2.FontSize = 10;
% cb2.Label.String = 'Degree Difference in Direction';
% cb2.Label.FontSize = 12;


%---------------------------------------Do the same thing for Average Wind Directions 
% ------------- Grab U AND V for Future and Historic 
temp_u_hist = nan_mean(u_avg(:,:,years1),3);
temp_v_hist = nan_mean(v_avg(:,:,years1),3);
wnddir_dif_matrixU = wrap2360(90-((180/pi)*atan2(temp_u_hist,temp_v_hist))+180); % Convert to Nautical 

temp_u_fut = nan_mean(u_avg(:,:,years2),3);
temp_v_fut = nan_mean(v_avg(:,:,years2),3);
wnddir_dif_matrixV = wrap2360(90-((180/pi)*atan2(temp_u_fut,temp_v_fut))+180); % Convert to Nautical 


% ------------- Plot U
row = 0;
col = 0;
a1 = axes('position',[p_left+col*(p_wid+p_spacing) p_bot+row*(p_height+p_spacing) p_wid p_height]);
pcolor(lon,lat,wnddir_dif_matrixU);
shading interp
hold on
plot(wa_lon, wa_lat,'k')
axis equal
xlim([min(lon(:)) -121.9])
ylim([min(lat(:)) 49.1])
caxis([200 340])
ylabel('Degrees Latitude')
xlabel('Degrees Longitude')

% ------------- Plot V
col = 1;
a2 = axes('position',[p_left+col*(p_wid+p_spacing) p_bot+row*(p_height+p_spacing) p_wid p_height]);
pcolor(lon,lat,wnddir_dif_matrixV);
shading interp
hold on
plot(wa_lon, wa_lat,'k')
axis equal
xlim([min(lon(:)) -121.9])
ylim([min(lat(:)) 49.1])
caxis([200 340])
xlabel('Degrees Longitude')
set(gca,'YTickLabel',[])

% ------------- Plot Difference
col = 2;
a3 = axes('position',[p_left+col*(p_wid+p_spacing) p_bot+row*(p_height+p_spacing) p_wid p_height]);
cmap = colormap(a3,mycolors);
dif_mat = wnddir_dif_matrixV - wnddir_dif_matrixU;
pcolor(lon,lat,dif_mat);
shading interp
hold on
plot(wa_lon, wa_lat,'k')
axis equal
xlim([min(lon(:)) -121.9])
ylim([min(lat(:)) 49.1])
xlabel('Degrees Longitude')
set(gca,'YTickLabel',[])
caxis([-10 10])

% ------------- Add Colorbars
cb1 = colorbar(a1,'Location','NorthOutSide','Position',[p_left+.01 .92 2*p_wid-.01 0.025]);
cb1.FontSize = 10;
cb1.Label.String = 'Wind Direction [degrees]';
cb1.Label.FontSize = 12;

cb2 = colorbar(a3,'Location','NorthOutSide','Position',[(p_left+2*p_wid+2*p_spacing)+.01 .92 p_wid-.02 0.025]);
cb2.FontSize = 10;
cb2.Label.String = 'Degree Difference in Direction';
cb2.Label.FontSize = 12;

% Save Figure
printFig(gcf,'Direction_DifferenceUV_Plot_MaxAndAvg',[14 8.5],'png',300)
