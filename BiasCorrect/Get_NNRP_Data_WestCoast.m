% This script goes through and grabs all the NNRP poitns within a window on
% the west coast mainly around WA for analysis 


clearvars
addpath C:\Functions_Matlab
model = 'NNRP';
switch model
    case 'NNRP'
        % Load Model domain grids
        cur_dir = pwd;
        fol_loc = 'E:\Abbas\Model_Met_Forcings\Wind_Products\ExtractNNRP\';
        cd(fol_loc)
        domaininfo = ncinfo('geo_em.d02.nc');
        lat_m = ncread('geo_em.d02.nc','XLAT_M');
        lon_m = ncread('geo_em.d02.nc','XLONG_M');
%         lat_u = ncread('geo_em.d02.nc','XLAT_U');
%         lon_u = ncread('geo_em.d02.nc','XLONG_U');
%         lat_v = ncread('geo_em.d02.nc','XLAT_V');
%         lon_v = ncread('geo_em.d02.nc','XLONG_V');
%         lat_c = ncread('geo_em.d02.nc','CLAT');     %Computational grid (same as mass grid)
%         lon_c = ncread('geo_em.d02.nc','CLONG');
        elev_m = ncread('geo_em.d02.nc','HGT_M');    %Height of topograpghy in meters MSL
        land_mask = ncread('geo_em.d02.nc','LANDMASK');
        land_mask = logical(land_mask);
%         [N, M] = size(lat_m);
        cd(cur_dir)
        
        % Rows and Columns to use for grabbing data of WA and surrounding
        % area
        r = 27:70;
        c = 50:112;
                
        %plot(lon_m(25:80,30:120),lat_m(25:80,30:120),'*')
        %bound = (lat_m >= 44 & lat_m <= 51 & lon_m >= -127 & lon_m <= -120)
        
        % Set times to load
        time_vec = datenum(1950,1,1):(1/4):datenum(2010,12,31,18,0,0);
        tic
        count = 0;
        
        slp = NaN(length(r),length(c),length(time_vec));
        wndspd = slp;
        wnddir = slp;
        surfp = slp;
        wind_u = slp;
        wind_v = slp;
        %sst = slp;
        %maxwndspd = slp;
        
        % Populate with DATA
        tic
        for tt=1:length(time_vec)
            % Load model prediction time step
            fol_loc = sprintf('E:/Abbas/Model_Met_Forcings/Wind_Products/NNRP/%d/',year(time_vec(tt)));
            file_name = sprintf('wrfoutp_d02_%04d%02d%02d%02d.nc',...
                year(time_vec(tt)),month(time_vec(tt)),day(time_vec(tt)),hour(time_vec(tt)));
            model_path = [fol_loc file_name];
            modelinfo = ncinfo(model_path);
            
            % Load in desired variables
            % time = ncread(model_path,'Time');
            surf_press = ncread(model_path,'PSFC');
            u10 = ncread(model_path,'U10');
            v10 = ncread(model_path,'V10');
%             maxwspd_mat = ncread(model_path,'MAXWSPD');
%             sst_mat = ncread(model_path,'SST')+273.15; %[C]
            Temp2m = ncread(model_path,'T2'); %[K]

            % Convert surface pressure to SLP
            surf_press_adj = convertToSLP(surf_press,elev_m,Temp2m);
            
            % Extract just neareast grid cells
            slp(:,:,tt) = double(surf_press_adj(r,c));
            surfp(:,:,tt) = double(surf_press(r,c));
            wind_u(:,:,tt) = double(u10(r,c));
            wind_v(:,:,tt) = double(v10(r,c));
            u10 = double(u10(r,c));
            v10 = double(v10(r,c));
%             sst(tt,:) = double(sst_mat(bound));
%             maxwspd(tt,:) = double(maxwspd_mat(bound));
            
            % Calc wind speed
            wndspd(:,:,tt) = sqrt(u10.^2 + v10.^2);
            wnddir_temp = (180/pi)*atan2(v10,u10);
            
            % Rotate winds to compass directions
            wnddir_temp = 90 - wnddir_temp;
            wnddir_temp(wnddir_temp<0)=wnddir_temp(wnddir_temp<0)+360;
            
            % switch to conventional coming from dir
            wnddir_temp = wnddir_temp+180;
            wnddir(:,:,tt) = wrap2360(wnddir_temp);
            
            if rem(tt,10000) == 0
                fprintf('Completed: %2.1f Percent\n',tt/length(time_vec))
            end
            clear time surf_press surf_press_adj Temp2m
        end
        toc
end
time = time_vec';

% Subsample
lat_m = lat_m(r,c);
lon_m = lon_m(r,c);
elev_m = elev_m(r,c);
land_mask = land_mask(r,c);
save('west_coast_NNRP.mat','wndspd','wnddir','wind_u','wind_v','time','slp','surfp','lat_m','lon_m','elev_m','land_mask')
toc