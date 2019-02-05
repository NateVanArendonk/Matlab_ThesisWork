% study the color content of brightest transect image data.

clear
load('OUTPUT/1444926217.Thu.Oct.15_16_23_37.GMT.2015.argus02b.cx.RMaxDay.mat')
i = 7;      % example sunny day image
cam = 1;

% extract an image, the R,G and B color planes. 
IInt = imread(FTPPath(RMaxDay.RMaxes(i).fnList(cam,:)));
Ig = double(rgb2gray(IInt));
I = double(IInt);

% pick a transect manually to span dune, sand, breakers.
f = ginput(2);          % pick a transect manually by two endpoints
U = (f(1,1)): f(2,1);
V = interp1([f(1,1) f(2,1)], [f(1,2) f(2,2)], U);
R = I(:,:,1);
G = I(:,:,2);
B = I(:,:,3);
Us = 1:size(R,2); Vs = 1: size(R,1);
IR = interp2(Us,Vs,R,U,V);      % grab R, G, B and gray shade intensities
IG = interp2(Us,Vs,G,U,V);
IB = interp2(Us,Vs,B,U,V);
Ig = interp2(Us,Vs,Ig,U,V);

% plot original image with line
figure(1); clf
imagesc(IInt); hold on
line([f(1,1) f(2,1)], [f(1,2) f(2,2)])
axis off; title(datestr(epoch2Matlab(RMaxDay.RMaxes(i).when)))

% who R, G and B channels
figure(2); clf
plot(IR, 'r'); hold on
plot(IG, 'g')
plot(IB, 'b')
xlabel('pixel along transect'); ylabel('intensity'); grid on
title(datestr(epoch2Matlab(RMaxDay.RMaxes(i).when)))

% show gray shade intensity and ratio of blue to red
figure(3); clf
plot(IB./IR)
hold on;
plot(Ig/max(Ig), 'k')
xlabel('pixel along transect'); ylabel('intensity')
legend('ratio blue/red', 'gray', 'location', 'southeast'); grid on
title(datestr(epoch2Matlab(RMaxDay.RMaxes(i).when)))

% create BR image directly, then sample
BR = B./R;
figure; imagesc(BR); caxis([0 2])

