clf 
clearvars
hindcastFolder = 'E:\Abbas\PS_COSMOS\Thesis_Modeling\LUT\Hindcast_Output\';

h = dir('E:\Abbas\PS_COSMOS\Thesis_Modeling\LUT\Hindcast_Output\*.mat');
tic
lx = [];
ly = [];
for ii = 1:length(h)
    load([hindcastFolder h(ii).name]);
    lx = vertcat(lx,lon);
    ly = vertcat(ly,lat);
    plot(lon,lat,'.')
    hold on 
    fprintf('Completed %d\n',ii)
end
toc
%% 


