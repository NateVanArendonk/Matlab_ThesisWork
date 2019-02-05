% paths to model directory
clearvars 
clc
model_folder = 'C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\Tacoma\Owen\XBeach\OB_Runs\RoughnessTesting\chezy\';
runs = dir('C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\Tacoma\Owen\XBeach\OB_Runs\RoughnessTesting\chezy\*'); % list of run folders 
runs(1:2) = [];



%% Will automate when the time comes 

run_xb = 'call "xbeach.exe"';

cd(path_model)
system(run_xb)

