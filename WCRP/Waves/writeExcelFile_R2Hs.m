clearvars
clc
F = dir('WCRP_R2_HsigOut/*.mat');
for ff = 1:length(F)
    W = load(['WCRP_R2_HsigOut/' F(ff).name]);
    fname = F(ff).name;
    Ind = strfind(fname,'.mat');
    fname = fname(1:Ind-1);
    fname = strcat(fname,'.xlsx');
%     fid = fopen(fname,'w');
%     fprintf(fid,'Hsig, Tp, R2_Low, R2_High,\n %1.2f, %1.2f, %1.2f, %1.2f',W.hs,W.tp,W.R2(1),W.R2(2));
%     fclose(fid)
    contents = {'Hsig','Tp','R2_Low','R2_High';W.hs,W.tp,W.R2(1),W.R2(2)};
    xlswrite(fname,contents)
    movefile(fname,'ExcelOut')
end