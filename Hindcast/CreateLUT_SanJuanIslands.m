% Develop LUT from Nested SWAN Runs
% - Code modified from Skagit LUT written by S.Crosby: N.V. 12/7/2017
% 1. Load  model run data
% 2. Determine grid cell or two from NNRP predictions to use
% 3. Create time series of wind speed, dir, and water level
% 4. Create interpolation from LUT given speed, dir and wl
% 5. Determine output points (or transects)
clearvars
addpath C:/Functions_Matlab


% SWAN runs location
fol_run = 'E:\Abbas\PS_COSMOS\Thesis_Modeling\SWAN\PS_RegionalModel\Model_Runs\NV_LUT\RES1\'; % Location of Swan Model Runs
kml_fol = 'E:\Abbas\PS_COSMOS\Thesis_Modeling\KML\PugetSoundShorline\SanJuanShorelines\'; % Location of shoreline KML to use
kml_list = dir('E:\Abbas\PS_COSMOS\Thesis_Modeling\KML\PugetSoundShorline\SanJuanShorelines\*.kml'); % List of all island alongshore transects
% kml_mask = 'C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\KML\Basin_Masks\LargeBasinMasks\South_PS_SWAN_Domain.kml'; % Used to determine which model to use

% Parameter Ranges
tide = [-2,-1,0,1,2,3,4,5.5];
speed = 0:5:30;
direc = -10:10:370;

% set minimum depth
min_depth = .1; %[m];


for kk = 1:length(kml_list)
    % Load alongshore offshore transect
    kml = kml_list(kk).name;
    extInd = strfind(kml,'.kml');
    tempName = kml(1:extInd-1);
    out_file = sprintf('%s_10m_LUT_offshore_extract.mat',tempName);
    if ~exist(['Offshore_LUT_Extracted\' out_file])
        fprintf('%s does not exist, making LUT', tempName)
        % Load KML Data - Basin Mask and Transect
        K = kml2struct([kml_fol kml]);
        % M = kml2struct(kml_mask);
        
        % Convert to utm
        [K.x,K.y] = deg2utm(K.Lat,K.Lon);
        
        % Create finely spaced alongshore transect to be derefined later
        x_e = [];
        y_e = [];
        for nn = 1:length(K.x)-1
            step = 1; % 1 meter
            [ line_x, line_y ] = createTransect( K.x(nn), K.y(nn), K.x(nn+1), K.y(nn+1), step );
            x_e = cat(1,x_e,line_x');
            y_e = cat(1,y_e,line_y');
        end
        % Get rid of any nans
        nan_inds = isnan(x_e);
        x_e(isnan(x_e)) = [];
        y_e(isnan(y_e)) = [];
        
        % Derefine to be specific resolution
        res = 50; %[m]
        ra = 0:res:length(x_e);
        ra(1) = 1;
        x_e = x_e(ra);
        y_e = y_e(ra);
        % Convert back to lat lon to find grid indices
        [lat,lon] = utm2deg(x_e,y_e,char(ones(length(x_e),1)*'10 T'));
        K.lon = lon;
        K.lat = lat;
        
        % Initialize by loading a single run
        tt = 1;
        ss = 3;
        dd = 10;
        fname = sprintf('SpatialJDF_s%02d_d%d_t%d.mat',speed(ss)*10,direc(dd),tide(tt)*10); % Loads in a low tide level
        data = load([fol_run fname]);
        i_x = zeros(length(K.lon),1);
        i_y = i_x;
        
        % Load in the depths for the jdf and ps model
        % ------------------------- JDF Model -------------------------------------
        dep_nm = 'E:\Abbas\PS_COSMOS\Thesis_Modeling\SWAN\PS_RegionalModel\INP_sph2\jdaf_new.BOT';
        grd_nm = 'E:\Abbas\PS_COSMOS\Thesis_Modeling\SWAN\PS_RegionalModel\INP_sph2\jdf_sph_swn.grd';
        
        % From .swn file
        mx = 659;
        my = 1195;
        % Load grid
        J = swan_io_grd('read',grd_nm,mx,my,99,99);
        
        % Load bottom
        S.mxinp = mx;
        S.myinp = my;
        S.idla = 4;
        S.nhedf = 0;
        S.fname1 = dep_nm;
        S.quantity = 'depth';
        J.Z = swan_io_bot('read',grd_nm,S);
        
        % Get Depths Specific Grid of interest
        data.Zp = J.Z';
        % Make a Non-Nan vector to find closest actual wave grid point to 
        nanInds = isnan(data.Hsig);
        tx = data.Xp;
        ty = data.Yp;
        tx(nanInds) = -126;
        ty(nanInds) = 48;
        
        % Find nearest grid cells for extraction
        for nn = 1:length(K.lon)
            [i_x(nn), i_y(nn)] = findNearestGridPoint(tx,ty,K.lon(nn),K.lat(nn));
        end
        
        % Size of grid
        [n, m] = size(data.Hsig);
        
        % Save output locations
        E.x = data.Xp(sub2ind([n m],i_x,i_y));
        E.y = data.Yp(sub2ind([n m],i_x,i_y));
        E.z = data.Depth(sub2ind([n m],i_x,i_y));
        
        % Check and make sure that we grabbed a deep enough point, if not, grab a
        % new point
        highEnd = 10;
        lowEnd = 40;
        for ii = 1:length(E.z)
            count = 0;
            % Original Values if this doesn't work below
            OGX = E.x(ii);
            OGY = E.y(ii);
            OGZ = E.z(ii);
            while E.z(ii) < highEnd || E.z(ii) > lowEnd
                depthGrid = zeros(3,3);
                xGrid = depthGrid;
                yGrid = depthGrid;
                ixGrid = depthGrid;
                iyGrid = depthGrid;
                % Order to work through wave grid
                xx = [-1 0 1 -1 0 1 -1 0 1];
                yy = [1 1 1 0 0 0 -1 -1 -1];
                % Grab depths at each point in box around current grid point
                % Note Depths are Positive!
                for jj = 1:3
                    for kk = 1:length(xx)
                        mxx = data.Xp(sub2ind([n, m],i_x(ii)+xx(kk),i_y(ii)+yy(kk)));
                        myy = data.Yp(sub2ind([n, m],i_x(ii)+xx(kk),i_y(ii)+yy(kk)));
                        mzz = data.Depth(sub2ind([n, m],i_x(ii)+xx(kk),i_y(ii)+yy(kk)));
                        xGrid(2+xx(kk),2+yy(kk)) = mxx;
                        yGrid(2+xx(kk),2+yy(kk)) = myy;
                        depthGrid(2+xx(kk),2+yy(kk)) = mzz;
                        ixGrid(2+xx(kk),2+yy(kk)) = i_x(ii)+xx(kk);
                        iyGrid(2+xx(kk),2+yy(kk)) = i_y(ii)+yy(kk);
                    end
                end
                % Grab the deepest depth point, and if it is less than 10, the loop
                % will stop, otherwise, it will keep on going
                % Lets first see if any are greater than or equal to the
                % criteria I set
                logMat = depthGrid >= highEnd;
                if sum(logMat(:)) >= 1
                    depthGrid(~logMat) = 999;
                    diffMat = ones(size(depthGrid))*highEnd;
                    diffMat = abs(depthGrid - diffMat);
                    [~,Id] = min(diffMat(:),[],'omitnan');
                    % Get Row and Column
                    [I_row, I_col] = ind2sub(size(depthGrid),Id);
                else % we need to relocate - find the closest point to my criteria
                    diffMat = ones(size(depthGrid))*highEnd;
                    diffMat = abs(depthGrid - diffMat);
                    [~,Id] = min(diffMat(:),[],'omitnan');
                    % Get Row and Column
                    [I_row, I_col] = ind2sub(size(depthGrid),Id);
                end
                % Change wave grid point to be new depth location
                E.x(ii) = xGrid(I_row,I_col);
                E.y(ii) = yGrid(I_row,I_col);
                E.z(ii) = depthGrid(I_row,I_col);
                i_x(ii) = ixGrid(I_row,I_col);
                i_y(ii) = iyGrid(I_row,I_col);
                count = count + 1;
                if count > 20
                    %                 if E.z(ii) > OGZ
                    %                     E.x(ii) = OGX;
                    %                     E.y(ii) = OGY;
                    %                     E.z(ii) = OGZ;
                    %                 end
                    if E.z(ii) > lowEnd
                        fprintf('WARNING - DEPTH IS GREATER THAN USER THRESHOLD - not an alarm just note it\n')
                    end
                    break
                end
            end
        end
        % See if there are any nan's
        tempHs = data.Hsig(sub2ind([n m],i_x,i_y))';
        nanInds = find(isnan(tempHs));
        if ~isempty(nanInds)
            fprintf('Error with %s - NaN points still exist in model extraction\n',tempName)
            saveNm = sprintf('%s_HAS_NAN.txt',tempName);
            save(saveNm);
        end
        % This should hopefully get rid of any NaN indices and be better in
        % the long run.
%         while ~isempty(nanInds)
%             for ii = 1:length(nanInds)
%                 count = 0;
%                 ind = nanInds(ii);
%                 % Original Values if this doesn't work below
%                 OGX = E.x(ind);
%                 OGY = E.y(ind);
%                 OGZ = E.z(ind);
%                 while E.z(ind) == OGZ || E.z(ind) < highEnd
%                     depthGrid = zeros(3,3);
%                     xGrid = depthGrid;
%                     yGrid = depthGrid;
%                     ixGrid = depthGrid;
%                     iyGrid = depthGrid;
%                     % Order to work through wave grid
%                     xx = [-1 0 1 -1 0 1 -1 0 1];
%                     yy = [1 1 1 0 0 0 -1 -1 -1];
%                     % Grab depths at each point in box around current grid point
%                     % Note Depths are Positive!
%                     for jj = 1:3
%                         for kk = 1:length(xx)
%                             mxx = data.Xp(sub2ind([n, m],i_x(ii)+xx(kk),i_y(ii)+yy(kk)));
%                             myy = data.Yp(sub2ind([n, m],i_x(ii)+xx(kk),i_y(ii)+yy(kk)));
%                             mzz = data.Depth(sub2ind([n, m],i_x(ii)+xx(kk),i_y(ii)+yy(kk)));
%                             xGrid(2+xx(kk),2+yy(kk)) = mxx;
%                             yGrid(2+xx(kk),2+yy(kk)) = myy;
%                             depthGrid(2+xx(kk),2+yy(kk)) = mzz;
%                             ixGrid(2+xx(kk),2+yy(kk)) = i_x(ii)+xx(kk);
%                             iyGrid(2+xx(kk),2+yy(kk)) = i_y(ii)+yy(kk);
%                         end
%                     end
%                     % Grab the deepest depth point, and if it is less than 10, the loop
%                     % will stop, otherwise, it will keep on going
%                     % Lets first see if any are greater than or equal to the
%                     % criteria I set
%                     logMat = depthGrid >= highEnd;
%                     if sum(logMat(:)) >= 1
%                         depthGrid(~logMat) = 999;
%                         diffMat = ones(size(depthGrid))*highEnd;
%                         diffMat = abs(depthGrid - diffMat);
%                         [~,Id] = min(diffMat(:),[],'omitnan');
%                         % Get Row and Column
%                         [I_row, I_col] = ind2sub(size(depthGrid),Id);
%                     else % we need to relocate - find the closest point to my criteria
%                         diffMat = ones(size(depthGrid))*highEnd;
%                         diffMat = abs(depthGrid - diffMat);
%                         [~,Id] = min(diffMat(:),[],'omitnan');
%                         % Get Row and Column
%                         [I_row, I_col] = ind2sub(size(depthGrid),Id);
%                     end
%                     % Change wave grid point to be new depth location
%                     E.x(ind) = xGrid(I_row,I_col);
%                     E.y(ind) = yGrid(I_row,I_col);
%                     E.z(ind) = depthGrid(I_row,I_col);
%                     i_x(ind) = ixGrid(I_row,I_col);
%                     i_y(ind) = iyGrid(I_row,I_col);
%                 end
%             end
%             count = count + 1;
%             % See if there are any nan's
%             tempHs = data.Hsig(sub2ind([n m],i_x,i_y))';
%             nanInds = find(isnan(tempHs));
%             if count > 20
%                 fprintf('Didnt Work, Quitting\n')
%                 break
%             end
%         end
        
        
        % Loop over parameterizations (takes 15 minutes)
        tic
        E.hs = NaN(length(i_x),length(tide),length(speed),length(direc));
        E.tp = E.hs;
        E.tm = E.hs;
        E.hdir = E.hs;
        
        for tt = 1:length(tide)
            for ss = 1:length(speed)
                for dd = 1:length(direc)
                    if speed(ss) ~= 0
                        fname = sprintf('SpatialJDF_s%02d_d%d_t%d.mat',speed(ss)*10,direc(dd),tide(tt)*10);
                        if direc(dd) == -10
                            ind = find(direc == 350);
                            fname = sprintf('SpatialJDF_s%02d_d%d_t%d.mat',speed(ss)*10,direc(ind),tide(tt)*10);
                        elseif direc(dd) == 360
                            ind = find(direc == 0);
                            fname = sprintf('SpatialJDF_s%02d_d%d_t%d.mat',speed(ss)*10,direc(ind),tide(tt)*10);
                        elseif direc(dd) == 370
                            ind = find(direc == 10);
                            fname = sprintf('SpatialJDF_s%02d_d%d_t%d.mat',speed(ss)*10,direc(ind),tide(tt)*10);
                        end
                        data = load([fol_run fname]);
                        dep = data.Depth(sub2ind([n m],i_x,i_y))';
                        dry_inds = dep < min_depth;
                        
                        % Get alongshore point wave info
                        E.hs(:,tt,ss,dd) = data.Hsig(sub2ind([n m],i_x,i_y))';
                        E.tp(:,tt,ss,dd) = data.RTpeak(sub2ind([n m],i_x,i_y))';
                        E.tm(:,tt,ss,dd) = data.Tm01(sub2ind([n m],i_x,i_y))';
                        E.hdir(:,tt,ss,dd) = data.Dir(sub2ind([n m],i_x,i_y))';
                        
                        % Set NaN dry indices
                        E.hs(dry_inds,tt,ss,dd) = 0;
                        E.tp(dry_inds,tt,ss,dd) = 0;
                        E.tm(dry_inds,tt,ss,dd) = 0;
                        E.hdir(dry_inds,tt,ss,dd) = 0;
                    else % Set values to zero if wind speed is less than 1
                        % Get alongshore point wave info
                        E.hs(:,tt,ss,dd) = 0;
                        E.tp(:,tt,ss,dd) = 0;
                        E.tm(:,tt,ss,dd) = 0;
                        E.hdir(:,tt,ss,dd) = 0;
                    end
                end
            end
            fprintf('Completed Tide level %d - Moving On\n',tide(tt))
        end
        toc
        
        E.tide = tide;
        E.speed = speed;
        E.direc = direc;
        E.lat = K.lat;
        E.lon = K.lon;
        E.x_u = x_e;
        E.y_u = y_e;
        E.gridX = data.Xp;
        E.gridY = data.Yp;
        
        save(out_file,'-struct','E','-v7.3')
        movefile(out_file,'Offshore_LUT_Extracted');
    end
    
end




