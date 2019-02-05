clearvars

addpath C:\Functions_Matlab
% runs location
fol_run = 'E:\Abbas\PS_COSMOS\Thesis_Modeling\SWAN\PS_RegionalModel\Model_Runs\NV_LUT\RES2\';

% Ranges for variables used
tide = [-2,-1,0,1,2,3,4,5.5];
speed = 5:5:30;
direc = 0:10:350;






%% Make a movie of figures - Doesn't work on Abbas, doesn't have software packages :(


v = VideoWriter('PS_Model_Spotlight','MPEG-4');
v.FrameRate = .8;
v.Quality = 75;
open(v)
figure
for tt = 1
    for ss = 3
        for dd = 1:length(direc)
            clf
            fname = sprintf('SpatialPS_s%d_d%d_t%d',speed(ss)*10,direc(dd),tide(tt)*10);
            data = load([fol_run fname]);
            pcolor(data.Xp,data.Yp,data.Hsig)
            shading flat
            caxis([0 2])
            cb = colorbar;
            cb.Label.String = 'Hsig [m]';
            xlabel('Degrees Longitude')
            ylabel('Degrees Latitude')
            axis equal
            frame = getframe(gcf);
            writeVideo(v,frame);
            
        end
    end
end
close(v)
