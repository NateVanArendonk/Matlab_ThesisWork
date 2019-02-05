clearvars
clc
F = dir('WCRP_WaveOut/*.mat'); % Load innames of all of the structures of circle wave points
for ff = 1:length(F)
    load(['WCRP_WaveOut/' F(ff).name]) % Load in each circle
    for hh= 1:length(H) % Go through and find the max yearly wave conditions from each level of the structure
        [M(hh).hs,I] = max(H(hh).hs_yr_avg);
        M(hh).tp = H(hh).tp_yr_avg(I);
        M(hh).x = H(hh).x(I);
        M(hh).y = H(hh).y(I);
        M(hh).z = H(hh).z(I);
        M(hh).depth = H(hh).depth(I);
    end
    % Destructify the structure variables
    hs = zeros(length(M),1);
    tp = hs;
    x = hs;
    y = hs;
    z = hs;
    depth = hs;
    % Populate with variables
    for ii = 1:length(M)
        hs(ii) = M(ii).hs;
        tp(ii) = M(ii).tp;
        x(ii) = M(ii).x;
        y(ii) = M(ii).y;
        z(ii) = M(ii).z;
        depth(ii) = M(ii).depth;
    end
    % Grab max from list of maxs and that will be max wave and location
    [~,I] = max(hs);
    hs = hs(I);
    x = x(I);
    y = y(I);
    z = z(I);
    depth = depth(I);
    
    % Perform linear regression to find associated Tp
    hsr = [];
    tpr = [];
    for ii = 1:length(H)
        hsr = vertcat(hsr,H(ii).hs_yr_avg);
        tpr = vertcat(tpr,H(ii).tp_yr_avg);
    end
    [p,~,~] = polyfit(hsr,tpr,1);
    % p(1) is the slope
    % p(2) is the intercept
    yfit = polyval(p,hsr); % Get y component of line
    yresid = tpr - yfit; % Compute Residuals
    SSresid = sum(yresid.^2); % Calculate residual sum of squares
    SStotal = (length(tpr)-1)*var(tpr); % Compute total sum of squares of y 
    rsq = 1 - SSresid/SStotal; % Calculate r-squared
    
    tp = p(1)*hs+p(1); % Solve linear equation for Tp given max Hsig
    saveNm = F(ff).name;
    ind = strfind(saveNm,'_HsOut');
    saveNm = saveNm(1:ind-1);
    saveNm = strcat(saveNm,'_HsMax.mat');
    save(saveNm,'hs','tp','x','y','z','depth','rsq');
    movefile(saveNm,'Wave_CircleYrMax')
end



