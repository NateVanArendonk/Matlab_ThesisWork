% Run LUT and Hindcast 
clearvars
clc
CreateLUT_PS;
movefile('SWAN_10m_PSLUT_offshore_extract.mat','Offshore_LUT_Extracted');
CreateLUT_JDF;
movefile('SWAN_10m_JDFLUT_offshore_extract.mat','Offshore_LUT_Extracted');
CreateLUT_SGA;
movefile('SWAN_10m_SGALUT_offshore_extract.mat','Offshore_LUT_Extracted');
CreateLUT_SanJuanIslands;
fprintf('San Juan Islands are Complete!\n')
CreateLUT_SouthPSIslands;
fprintf('South Puget Sound Islands are Complete!\n')
CreateHindcast;
GetHsig;
return
%% old transfering code 

if exist('WaldronIsland_10m_LUT_offshore_extract.mat')
    movefile('WaldronIsland_10m_LUT_offshore_extract.mat','Offshore_LUT_Extracted');
end
if exist('Cypress_Island_10m_LUT_offshore_extract.mat')
    movefile('Cypress_Island_10m_LUT_offshore_extract.mat','Offshore_LUT_Extracted');
end
if exist('LopezIsland_10m_LUT_offshore_extract.mat')
    movefile('LopezIsland_10m_LUT_offshore_extract.mat','Offshore_LUT_Extracted');
end
if exist('OrcasIsland_10m_LUT_offshore_extract.mat')
    movefile('OrcasIsland_10m_LUT_offshore_extract.mat','Offshore_LUT_Extracted');
end
if exist('SanJuanIsland_10m_LUT_offshore_extract.mat')
    movefile('SanJuanIsland_10m_LUT_offshore_extract.mat','Offshore_LUT_Extracted');
end
if exist('StuartIsland_10m_LUT_offshore_extract.mat')
    movefile('StuartIsland_10m_LUT_offshore_extract.mat','Offshore_LUT_Extracted');
end
if exist('SuciaIsland_10m_LUT_offshore_extract.mat')
    movefile('SuciaIsland_10m_LUT_offshore_extract.mat','Offshore_LUT_Extracted');
end
if exist('Sinclair_Island_10m_LUT_offshore_extract.mat')
    movefile('Sinclair_Island_10m_LUT_offshore_extract.mat','Offshore_LUT_Extracted');
end
if exist('Guemes_Island_10m_LUT_offshore_extract.mat')
    movefile('Guemes_Island_10m_LUT_offshore_extract.mat','Offshore_LUT_Extracted');
end
if exist('Anderson_Island_10m_LUT_offshore_extract.mat')
    movefile('Anderson_Island_10m_LUT_offshore_extract.mat','Offshore_LUT_Extracted');
end
if exist('FoxIsland_10m_LUT_offshore_extract.mat')
    movefile('FoxIsland_10m_LUT_offshore_extract.mat','Offshore_LUT_Extracted');
end
if exist('HarstineIsland_10m_LUT_offshore_extract.mat')
    movefile('HarstineIsland_10m_LUT_offshore_extract.mat','Offshore_LUT_Extracted');
end
if exist('Vashon_10m_LUT_offshore_extract.mat')
    movefile('Vashon_10m_LUT_offshore_extract.mat','Offshore_LUT_Extracted');
end