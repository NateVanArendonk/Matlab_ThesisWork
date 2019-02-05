%% Load in variables of interest
addpath C:\Functions_Matlab
clear all
close all
clc

path = 'C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\Tacoma\Owen\XBeach\OwenBeach_R2\model_runs\OB_13_SLR\';
% path = 'C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\Tacoma\Owen\XBeach\OwenBeach_R2\model_runs\OB_Contemporary\';
grid_path = 'C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\Tacoma\Owen\XBeach\OwenBeach_R2\XbeachGrids\';
filename = 'xboutput.nc';

xb_info = ncinfo([path filename]);

% Load Grid Info
X.x = load([grid_path 'OB1_Xbeach_x.grd']);
X.y = load([grid_path 'OB1_Xbeach_y.grd']);
X.z = load([grid_path 'OB1_Xbeach_z.grd']);

X.point_zs = squeeze(ncread([path filename],'point_zs')); % Runup
X.point_xz = squeeze(ncread([path filename],'point_xz'));
X.point_yz = squeeze(ncread([path filename],'point_yz'));
X.zs = squeeze(ncread([path filename],'zs')); % Water Level

X.Hs = 4*sqrt(var(X.zs(:,1:end-1)'));
rho = 1027; %[kg/m3]
X.E = ((X.Hs.^2)*rho*9.81)/8; % Correct?

% Make cross shore 'S' transect
X.s = sqrt(X.x.^2+X.x.^2);
X.s = X.s - min(X.s);

% Make Cross shore 'S' transect for point Runup gauge data - BROKEN
X.point_s = sqrt(X.point_xz.^2+X.point_yz.^2);
X.point_s = X.point_s - min(X.point_s);
%%
v = VideoWriter('OwenBeach_1m4secWaves','MPEG-4');
v.FrameRate = 5;
v.Quality = 75;
open(v)
figure
for ii = 1:size(X.zs,1)
    clf
    fill([X.s X.s(end) X.s(1)], [X.zs(:,ii)' X.z(1) X.z(1)],'b','LineStyle','none') % Plot polygon of water
    hold on
    fill([X.s X.s(end) X.s(end)], [X.z X.z(end) X.z(1)],[1 .9 .4],'LineStyle','none') % Plot polygon of land
    grid on 
    hold off
    title(['t = ' sprintf('%04d',ii) ' seconds'])
    xlabel('Cross Shore [m]')
    ylabel('Wave height [m]')
    axis([min(X.s) max(X.s) -2.3 6])
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v);
