%runup_calculations_TAW_snsh


%% Barrier/Cliff Runup Reduction Factors - TAW Approach %%
%
%%% Roughness reduction factor (Eqn D.4.5.20) %%%
% Roughness factor (0.9) for vegetated (grass ~3 cm) slope. 
% See van der Meer and Janssen (1994) and CEM (2003)
% Roughness factor (0.6) for thickly vegetated slope 
% (i.e. salal/gorse etc). Dr. Bill McDougal, Coastal Engineer 
% (personal communication, April 2010)
% yr = 0.6;                                                   
% Clean up roughness factor
yr(isnan(yr)) = 1;                       % v0.1.1
clear yb                                    % model crashes if not
%%% Berm width reduction factor calculation (Eqn D.4.5-21) %v0.1, v0.1.5

if isnan(yb_or(aa))
    
    clear Lberm
    
    %Alternative LBerm calc (based on dh+/-Hmo) - added NTC 12 Sep 2013
    tmptop = Tide+Hmo; %added NTC 26 Sep 2013
    tmpbot = Tide-Hmo;

    
    for ii = 1:length(Hmo)

        if isnan(Hmo(ii)) == 0
            
            % Xtop (Location of STK_TWL)
            xint = curveintersect(x,y,...
                [x(1) x(end)+10],[tmptop(ii) tmptop(ii)]);
            Xtop = min(xint);

            % Xbottom (Location of Tide + STK_setup)
            xint = curveintersect(x,y,...
                [x(1) x(end)+10],[tmpbot(ii) tmpbot(ii)]);
            Xbottom = max(xint);

            Lberm(ii) = Xbottom - Xtop;
            
        elseif isnan(Hmo(ii)) == 1 %added NTC 13 Sep 2013
            
            Lberm(ii) = NaN;
        end
    
    end
    
% Determine xberm (Equation 4.5-21 & 13 (TAW report p17))JA v0.1.6

        if (1.1*(setup + STK_swash)+Tide)>-dh(aa)>0
            xberm=1.1*(setup + STK_swash)+Tide;
        elseif 2*Hmo>dh(aa)>=0
            xberm=2*Hmo;
        elseif or(-dh(aa)>=(1.1*(setup + STK_swash)+Tide), dh(aa)>=2*Hmo)
            xberm=1;
        end

    % Compute the runup reduction factor
    yb = 1-(Bw(aa)./(2.*Lberm))'.*(1+cos((pi.*dh(aa))./xberm));
    yb(isnan(yb)) = 1;        % No reduction factor otherwise
else
    % If the transect is not identified as a Berm then use the table value.
    yb(1:length(STK_TWL),1) = yb_or(aa);
end

% According to the guidelines the yb must be between 0.6 and 1.
yb(yb<0.6) = 0.6;
yb(yb>1) = 1;
    
%%% Wave direction(Eqn D.4.5-38)

    Beta = abs(WDir - SA(aa));               % v0.1.1 changed SA to SA(aa)
    yB = 1 - 0.0033*Beta;
    yB(Beta>=80) = 1 - 0.0033*80;
%% TAW Approach With Snsh %%  
% v0.1.3

% Compute the local slope FIRST ITERATION
% The slope will be computed between the TWL computed with stockdon and the
% level of the static setup + tide. 
Ytop2 = Tide + (1.5*Hmo);
% If profile is overtopped then pick the top of the dune (minus 1cm for
% interpolation purposes, see curveintersect)
Ytop2(Ytop2 > Dc(aa)) = Dc(aa) - 0.01;
% The bottom part of the slope is given by
Ybottom2 = Tide - (1.5*Hmo);

% Preallocate x variables
Xtop2 = nan(size(STK_TWL));
Xbottom2 = nan(size(STK_TWL));

% Find the location of Xtop and Xbottom
for ii = 1:length(STK_TWL)

    if flag_Hmo(ii) == 0 && isnan(Hmo(ii)) == 0
        % Xtop (Location of STK_TWL)
        xint = curveintersect(x,y,...
            [x(1) x(end)+10],[Ytop2(ii) Ytop2(ii)]);
        Xtop2(ii) = max(xint);

        % Xbottom (Location of Tide + STK_setup)
        xint = curveintersect(x,y,...
            [x(1) x(end)+10],[Ybottom2(ii) Ybottom2(ii)]);
        Xbottom2(ii) = max(xint);
    end
        
end

% Compute slope where the swash acts
S_TAW_Snsh = (Ytop2 - Ybottom2)./abs(Xtop2-Xbottom2);

%%%%%% RUNUP EQUATION - DETERMINISTIC APPROACH (includes +1STDEV factor of safety) %%%%%%
% Eqn 5.4 - Eurotop (2007): yr=roughness, yb=berm factor, yB = wave direction
Ib = S_TAW_Snsh./(sqrt(Hmo./L));                                % Irribarren number based on 2ND slope estimate - Eqn D.4.5-8
ybIb = yb.*Ib;                                             
TAW_R_Snsh = zeros(1,length(ybIb))';
for j = 1:length(ybIb)
    if ybIb(j)>=0 && ybIb(j)<1.8             
        TAW_R_Snsh(j) = Hmo(j)*1.75*yr(aa)*yb(j).*yB(j).*Ib(j);                      
    elseif ybIb(j)>=1.8
        TAW_R_Snsh(j) = Hmo(j)*yr(aa).*yB(j).*(4.3-(1.6./sqrt(Ib(j))));
    end
end

% COMBINE RUNUP VALUES TO CALCULATE COMBINED TWL %
R_final_Snsh = (TAW_R_Snsh + 1.1*STK_setup);    % Runup  (v0.1.2)
TAW_TWL_Snsh = (R_final_Snsh+Tide);                 % Total Water Level

% Compute the local slope SECOND ITERATION ~~~~~~~~~
% The slope will be computed between the TWL computed with stockdon and the
% level of the static setup + tide. 
Ytop3 = TAW_TWL_Snsh;
S_TAW_Snsh_flag(1:length(STK_TWL),1) = 1;
% If profile is overtopped then pick the top of the dune (minus 1cm for
% interpolation purposes, see curveintersect)
Ytop3(Ytop3 > Dc(aa)) = Dc(aa) - 0.01;
S_TAW_Snsh_flag(Ytop3 > Dc(aa)) = 2;

% Preallocate x variables
Xtop3 = nan(size(STK_TWL));

% Find the location of Xtop and Xbottom
for ii = 1:length(STK_TWL)

    if flag_Hmo(ii) == 0 && isnan(Hmo(ii)) == 0
        % Xtop (Location of STK_TWL)
        xint = curveintersect(x,y,...
            [x(1) x(end)+10],[Ytop3(ii) Ytop3(ii)]);
        Xtop3(ii) = max(xint);
    end
        
end

% Compute slope where the swash acts
S_TAW_Snsh_2 = (Ytop3 - Ybottom2)./abs(Xtop3-Xbottom2);


%%%%%% RUNUP EQUATION - DETERMINISTIC APPROACH (includes +1STDEV factor of safety) %%%%%%
% Eqn 5.4 - Eurotop (2007): yr=roughness, yb=berm factor, yB = wave direction
Ib = S_TAW_Snsh_2./(sqrt(Hmo./L));                                % Irribarren number based on 2ND slope estimate - Eqn D.4.5-8
ybIb = yb.*Ib;                                             
TAW_R_Snsh_2 = zeros(1,length(ybIb))';
for j = 1:length(ybIb)
    if ybIb(j)>=0 && ybIb(j)<1.8             
        TAW_R_Snsh_2(j) = Hmo(j)*1.75*yr(aa)*yb(j).*yB(j).*Ib(j);                      
    elseif ybIb(j)>=1.8
        TAW_R_Snsh_2(j) = Hmo(j)*yr(aa).*yB(j).*(4.3-(1.6./sqrt(Ib(j))));
    end
end

% COMBINE RUNUP VALUES TO CALCULATE COMBINED TWL %

R_final_Snsh = (TAW_R_Snsh_2 + 1.1*STK_setup);   % Runup  (v0.1.2)
TAW_TWL_Snsh = (R_final_Snsh + Tide);                % Total Water Level
