clearvars
load('MaxTWL_RustonWay.mat');
D = load('RustonSubset_Gridded_Coned.mat');

swl = E(189).wl;
swlS1 = swl + (1.3*0.3048);
swlS2 = swl + (2.8*0.3048);
clear E
return
%%
c = contour(D.X,D.Y,D.Z,[swlS2 swlS2]);
c(:,1) = []; 
x1 = c(1,:);
y1 = c(2,:);
bID = find(x1 < 1000);
x1 = x1(1:bID-1);
y1 = y1(1:bID-1);

%%

save('OB_RW_SWL_37231m_SLR280','x1','y1')