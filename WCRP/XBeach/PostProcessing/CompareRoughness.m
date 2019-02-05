clearvars
clc

r1 = load('OB1_R2_rough75.mat');
r2 = load('OB1_R2_rough85.mat');
r3 = load('OB1_R2_rough95.mat');

% ----------- Load Waves
W = load('C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\Hindcast\Ruston_LUT\OwenRuston_wave_hindcast.mat');

% ----------- Load DEM
E = load('E:\Abbas\Modeling Resources\PS_DEM\Ruston_Way\OwenBeach_CONED_Subset.mat');

% ----------- Make a gridded interped surface of DEM
myX = E.x_s(1):1:E.x_s(end);
myY = E.y_s(end):1:E.y_s(1);
tic
[mx,my] = meshgrid(myX,myY);
mz = griddata(E.x_s,E.y_s,E.z_s,mx,my);
toc

% ---------- Load FEMA 100 yr flood
F = load('C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\Tacoma\Ruston_1Pct_Annual_Chance_Flood_Line');
% convert to UTM
[F.x,F.y] = deg2utm(F.lat,F.lon);

% ---------- Load KML Transect 
kml_fol = 'C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\KML\RunupTransects\OwenBeach\';
kml_nm = 'OB1.kml';
temp = kml2struct([kml_fol kml_nm]); % load in KML of Transect
T.lat = temp.Lat;
T.lon = temp.Lon;
[T.x_utm,T.y_utm] = deg2utm(T.lat,T.lon); % Convert to utm
T.name = kml_nm;
clear kml_fol kmls kml_nm

% Find closest Fema point to end of transect 
dist_list = sqrt((T.x_utm(2) - F.x).^2 + (T.y_utm(2)-F.y).^2);
[~,I] = min(dist_list);
% Grab elevation at that point
[r,c] = findNearestGridPoint(mx,my,F.x(I),F.y(I));
F100 = mz(r,c);

%% % Grab max and TWL

t1 = r1.twl+r1.V2;
t2 = r2.twl+r2.V2;
t3 = r3.twl+r3.V2;

[m1,I1] = max(t1);
[m2,I2] = max(t2);
[m3,I3] = max(t3);


%% Plot Timeseries of R2 around Max R2 
clf
st = I1-50;
ed = I1+50;
b = parula(10);
l1 = plot(W.time(st:ed),r1.V2(st:ed),'o','Color',b(1,:),'MarkerFaceColor',b(1,:));
hold on 
l2 = plot(W.time(st:ed),r2.V2(st:ed),'o','Color',b(4,:),'MarkerFaceColor',b(4,:));
l3 = plot(W.time(st:ed),r3.V2(st:ed),'o','Color',b(7,:),'MarkerFaceColor',b(7,:));

plot(W.time(st:ed),r1.V2(st:ed),'Color',b(1,:))
hold on 
plot(W.time(st:ed),r2.V2(st:ed),'Color',b(4,:))
plot(W.time(st:ed),r3.V2(st:ed),'Color',b(7,:))

grid on 
lgd = legend([l1,l2,l3],'Roughness: 0.75','Roughness: 0.85', 'Roughness: 0.95','Location','Northeast');
lgd.FontSize = 14;

avg78 = mean(r2.V2 - r1.V2); avg78 = sprintf('Avg. Increase .75 to .85 = %.2fm',avg78);
avg89 = mean(r3.V2 - r2.V2); avg89 = sprintf('Avg. Increase .85 to .95 = %.2fm',avg89);
max78 = max(r2.V2 - r1.V2); max78 = sprintf('Max Increase .75 to .85 = %.2fm',max78);
max89 = max(r3.V2 - r2.V2); max89 = sprintf('Max Increase .85 to .95 = %.2fm',max89);

text(W.time(st+58),0.69,avg78,'FontSize',14)
text(W.time(st+58),0.65,avg89,'FontSize',14)
text(W.time(st+4),0.69,max78,'FontSize',14)
text(W.time(st+4),0.65,max89,'FontSize',14)

xlabel('Time [Date Number]')
ylabel('R2 Height [m]')
set(gca,'FontSize',14)
printFig(gcf,'Runup_Roughness_Timeseries',[15,11],'png',600)

%% Plot of varying roughness spatial map 
b = parula(10);
% ---------------- Roughness .75
r75 = [max(t1) max(t1)];
cr75 = contour(mx,my,mz,r75,'b');
cr75(:,1) = [];
temp = cr75; clear cr75
cr75.x = temp(1,:);
cr75.y = temp(2,:);

% ------------ Roughness .85
r85 = [max(t2),max(t2)];
cr85 = contour(mx,my,mz,r85,'m');
cr85(:,1) = [];
temp = cr85; clear cr85
cr85.x = temp(1,:);
cr85.y = temp(2,:);

% ------------ Roughness .95
r95 = [max(t3) max(t3)];
cr95 = contour(mx,my,mz,r95,'m');
cr95(:,1) = [];
temp = cr95; clear cr95
cr95.x = temp(1,:);
cr95.y = temp(2,:);


clf
 % Load in Geotiff
IM = load('C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\GeoTiffs\OwenBeach\OB_utm_zoom2.mat');
imagesc(IM.xm,IM.ym,IM.im)
set(gca,'ydir','normal')
hold on

% Add Contours 
contour_levels = (0:.25:8);
[C,h] = contour(mx,my,mz,contour_levels); % use contourf to plot filled in bathy 
% c = colorbar;
% c.Label.String = 'Elevation [m - NAVD88]';
h.LineColor = [.7 .7 .7];
xlabel('Northing [m]');
ylabel('Easting [m]');

% Label Contours 
v = 0:1:8;
clabel(C,h,v)

% Add Transect 
plot(T.x_utm,T.y_utm,'w')

% Add Fema 
f = plot(F.x,F.y,'r','LineWidth',1.5);

% Add Runup values 
cr75 = plot(cr75.x(1:end-100),cr75.y(1:end-100),'Color',b(1,:),'LineWidth',1.5);
cr85 = plot(cr85.x(1:end-100),cr85.y(1:end-100),'Color',b(4,:),'LineWidth',1.5);
cr95 = plot(cr95.x(1:end-100),cr95.y(1:end-100),'Color',b(7,:),'LineWidth',1.5);


xlim([5.3560*10^5, 5.357*10^5])
ylim([5.24001*10^6, 5.24008*10^6])

set(gca,'FontSize',14)

lgd = legend([cr75,cr85,cr95,f],'Roughness: .75','Roughness: .85','Roughness: .95','FEMA 100yr Flood','Location','NorthEast');
lgd.FontSize = 14;

% Add Text
xx = 5.3568*10^5;
yy = 5.24006*10^6;
txt = 'Along Transect:';
text(xx,yy,txt,'FontSize',14,'Color','w')
xx = 5.3568*10^5;
yy = 5.240058*10^6;
txt = sprintf('FEMA: %.3fm',F100);
text(xx,yy,txt,'FontSize',14,'Color','w')
xx = 5.3568*10^5;
yy = 5.240056*10^6;
txt = sprintf('R75: %.3fm',max(t1));
text(xx,yy,txt,'FontSize',14,'Color','w')
xx = 5.3568*10^5;
yy = 5.240054*10^6;
txt = sprintf('R85: %.3fm',max(t2));
text(xx,yy,txt,'FontSize',14,'Color','w')
xx = 5.3568*10^5;
yy = 5.240052*10^6;
txt = sprintf('R95: %.3fm',max(t3));
text(xx,yy,txt,'FontSize',14,'Color','w')

printFig(gcf,'FEMAvsVaryingRoughness',[15 11],'png',300)