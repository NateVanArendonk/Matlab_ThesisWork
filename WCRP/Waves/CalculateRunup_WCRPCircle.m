clear
clc
F = dir('Wave_CircleYrMax/*.mat');
for ff = 1:length(F)
    W = load(['Wave_CircleYrMax/' F(ff).name]);
    slopes = [.1 .3];
    roughness = [.7 1]; % Per TAW, 0.7 is single layer of armor rock 
    R2 = zeros(length(slopes),length(roughness));
    for ss = 1:length(slopes)
        for rr = 1:length(roughness)
        R2(ss,rr) = VanderMeer2002_WCRPAdapted_R2(W.hs,W.tp,slopes(ss),W.depth,-1,0,roughness(rr));    
        end
    end
    R2 = [min(min(R2)) max(max(R2))];
    hs = W.hs;
    tp = W.tp;
    x = W.x;
    y = W.y;
    z = W.z;
    depth = W.depth;
    rsq = W.rsq;
    % Save  
    saveNm = F(ff).name;
    ind = strfind(saveNm,'_HsMax');
    saveNm = saveNm(1:ind-1);
    saveNm = strcat(saveNm,'_R2HsOut.mat');
    save(saveNm,'hs','tp','x','y','z','depth','rsq','R2');
    movefile(saveNm,'WCRP_R2_HsigOut')
end
