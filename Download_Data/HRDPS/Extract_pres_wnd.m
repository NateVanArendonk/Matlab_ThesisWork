% Extract pressure and winds from Candadian Weather Predictions
%   Input time, forecast hour, and locations to extract
%   Required Files
%       'lamwestpoints' - contains lat/lon positions of grid cells
%       'CMC_hrdps_west_XXX_MSL_0_ps2.5km_...' - prediction files
%           And prediction files locations
%       'etopo1_bedrock_PugetSound.nc' - Contains EPOTO1 bedrock

% Candadian Weather Data
%   https://weather.gc.ca/grib/grib2_HRDPS_HR_e.html
%   Using the LAMWEST polar-stereographic grid
%   Grid point lat/lon values in ASCII file
% Description:
%     ni	685 (x)
%     nj	485 (y)
%     resolution at 60° N	2.5 km
%     coordinate of first grid point	44.6922° N  129.9011° W
%     (i,j) coordinate of North Pole	(578.0, 1900.0)
%     grid orientation (with respect to j axis)	-113.0°
% Geographical coordinates for grid points saved in "lamwestpoints"
%     Formatted as i, j, lat, lon
% Missing Files notes
%  ..._2017020300_P000-00.grib2 is missing and strangely 2017020218 is in
%  the folder you expect to find it

clearvars
echo OFF ALL

grid_rotation = -113;

%%%%% Output Variables %%%%%%
%prtmsl - Pressure_reduced_to_MSL_msl [Pa]
%wind_speed - [m/s]
%wind_dir - [deg] Compass, degrees easting from North
%   Variable Array Format: (time, forecast hour, location)

%%%%% Input Variables %%%%%%

% Define sites of interest to extract
fname = 'D:\GoogleDrive\CoastalChangeModeling\Workplan\WeatherStationMeta.csv';
Station = readStationMeta( fname );

% Locations to extract (can be array)
input_loc_lat = Station.lat;
input_loc_lon = Station.lon;
input_loc_name = Station.shortName;

% Remove any locations without lat
I = find(isnan(input_loc_lat));
input_loc_lat(I) = [];
input_loc_lon(I) = [];
input_loc_name(I) = [];

% Grib data folder location
fol_data = '../CANADA/';

% Output Folder Location
fol_output = 'Output';

% Model Grid Domain Size (N x M), hardcoded for now
N = 485;
M = 685;

%%%%%%%%%%%%%%%%%%%%%% INITIALIZE %%%%%%%%%%%%%%%%%%%%

% load lat/lon locations for x and y coordinates
fid = fopen('lamwestpoints');
temp = textscan(fid,'%d %d %f %f');
fclose(fid);
grid_i = temp{1};
grid_j = temp{2};
grid_lat = temp{3};
grid_lon = temp{4};

% Load rotations to earth-relative
%   From Bert Rubash, derived from lat/lon using python 
rotation = textread('Rotations.dat');


% Convert list of lat/lon to grid same as paramters
grid_LAT = reshape(grid_lat,[M N])';
grid_LON = reshape(grid_lon,[M N])';

%%%%%%%%%%%% Here loop over time and forecast steps %%%%%%%%%%%%%%%%
time_loop = datenum(2016,02,02):(1/4):datenum(2017,03,22);
forecast_loop = 0:47;

% Find nearest grid cell locations to each input_loc
for nn = 1:length(input_loc_lat)
    min_dist_to_gridcell = 3000; %meters;
    [ i_vec(nn,:), j_vec(nn,:) ] = findNearestGridPos( grid_lat, grid_lon, grid_i, grid_j, input_loc_lat(nn), input_loc_lon(nn), min_dist_to_gridcell );
end

%%%%%%%%%%%%%%% UNCOMMENT FOR Test Loop %%%%%%%%%%%%%%%%%%%%
%time_loop = [datenum(2016,03,01) datenum(2016,03,02)];
%forecast_loop = 0:47;
%%%%%%%%%%%%%%% UNCOMMENT FOR Test Loop %%%%%%%%%%%%%%%%%%%%

% Preallocate
output_prtmsl = zeros(length(time_loop),length(forecast_loop),length(input_loc_lat));
output_wind_speed = output_prtmsl;
output_wind_dir = output_prtmsl;
output_prtmsl_four = zeros(length(time_loop),length(forecast_loop),length(input_loc_lat),4);
output_wind_speed_four = output_prtmsl_four;
output_wind_dir_four = output_prtmsl_four;

for tt = 1:length(time_loop)
    tic
    % Model runtime (every 6 hours, hr = 00, 06, 12 18)
    input_yr = year(time_loop(tt));
    input_mo = month(time_loop(tt));
    input_day = day(time_loop(tt));
    input_hr = hour(time_loop(tt));
    
    % Forecast hour (0, 1, 2 ... 47)
    for ff=1:length(forecast_loop)
        input_forecast_hour = forecast_loop(ff);
        
        % Set time stamp
        time = datenum(input_yr,input_mo,input_day,input_hr,0,0);
        
        % Folder Name
        fol_name = sprintf('downloaded_canadian_grib_files_%04d%02d%02d/',input_yr,input_mo,input_day);
        
        % Get pressure reduced to msl
        fname = sprintf('CMC_hrdps_west_PRMSL_MSL_0_ps2.5km_%04d%02d%02d%02d_P%03d-00.grib2',...
            input_yr,input_mo,input_day,input_hr,input_forecast_hour);
        [ prtmsl ] = retrieveGriddedData( fname, fol_name, fol_data );
        
        if ~isnan(prtmsl)
            
            % Get wind u
            fname = sprintf('CMC_hrdps_west_UGRD_TGL_10_ps2.5km_%04d%02d%02d%02d_P%03d-00.grib2',...
                input_yr,input_mo,input_day,input_hr,input_forecast_hour);
            [ wind_u ] = retrieveGriddedData( fname, fol_name, fol_data );
            wind_u = reshape(wind_u,[N M]);
            
            % Get wind v
            fname = sprintf('CMC_hrdps_west_VGRD_TGL_10_ps2.5km_%04d%02d%02d%02d_P%03d-00.grib2',...
                input_yr,input_mo,input_day,input_hr,input_forecast_hour);
            [ wind_v ] = retrieveGriddedData( fname, fol_name, fol_data );
            wind_v = reshape(wind_v,[N M]);
                                    
            % Calc wind speed
            wind_speed = sqrt(wind_u.^2 + wind_v.^2);
            wind_dir = (180/pi)*atan2(wind_v,wind_u);
            
            % Rotate winds to compass directions
            wind_dir = 90 - wind_dir;
            wind_dir(wind_dir<0)=wind_dir(wind_dir<0)+360;
            
            % switch to conventional coming from dir
            wind_dir = wind_dir+180; 
            wind_dir = wrapTo360(wind_dir);
                        
            % Interpolate to locations where data desired (SLOW!!)
%             debug.interp_prtmsl = griddata(grid_LON,grid_LAT,prtmsl,input_loc_lon,input_loc_lat);
%             debug.interp_wind_speed = griddata(grid_LON,grid_LAT,wind_speed,input_loc_lon,input_loc_lat);
%             debug.interp_wind_dir = griddata(grid_LON,grid_LAT,wind_dir,input_loc_lon,input_loc_lat);
            
            % Extract weather for each location
            % Use closest cell           
            interp_prtmsl = prtmsl(sub2ind([N M],j_vec(:,1),i_vec(:,1)));
            interp_wind_speed = wind_speed(sub2ind([N M],j_vec(:,1),i_vec(:,1)));
            
            % For wind_Dir Rotate to earth-relative coordinates
            interp_wind_dir = wind_dir(sub2ind([N M],j_vec(:,1),i_vec(:,1)))...
                - rotation(sub2ind([N M],j_vec(:,1),i_vec(:,1)));
               
%           % Save nearest 4 cells
            for nn=1:4
                four_prtmsl(:,nn) = prtmsl(sub2ind([N M],j_vec(:,nn),i_vec(:,nn)));
                four_wind_speed(:,nn) = wind_speed(sub2ind([N M],j_vec(:,nn),i_vec(:,nn)));
                four_wind_dir(:,nn) = wind_dir(sub2ind([N M],j_vec(:,nn),i_vec(:,nn)))...
                - rotation(sub2ind([N M],j_vec(:,nn),i_vec(:,nn)));
            end
            
            % Reshape and save in time and forecast hour
            output_prtmsl(tt,ff,:) = permute(interp_prtmsl,[3 2 1]);
            output_wind_speed(tt,ff,:) = permute(interp_wind_speed,[3 2 1]);
            output_wind_dir(tt,ff,:) = permute(interp_wind_dir,[3 2 1]);
            
            % Reshape and save in time and forecast hour (4 cells)
            output_prtmsl_four(tt,ff,:,:) = permute(four_prtmsl,[4 3 1 2]);
            output_wind_speed_four(tt,ff,:,:) = permute(four_wind_speed,[4 3 1 2]);
            output_wind_dir_four(tt,ff,:,:) = permute(four_wind_dir,[4 3 1 2]);
            
        else
            % If no file available fill with NaN
            output_prtmsl(tt,ff,:) = NaN(1,1,length(input_loc_lat));
            output_wind_speed(tt,ff,:) = NaN(1,1,length(input_loc_lat));
            output_wind_dir(tt,ff,:) = NaN(1,1,length(input_loc_lat));
            output_prtmsl_four(tt,ff,:,:) = NaN(1,1,length(input_loc_lat),4);
            output_wind_speed_four(tt,ff,:,:) = NaN(1,1,length(input_loc_lat),4);
            output_wind_dir_four(tt,ff,:,:) = NaN(1,1,length(input_loc_lat),4);
        end
    end
    fprintf(1,'Completed %s\n',datestr(time));
    toc
end


%%
% Format and save variables for saving

for nn = 1:length(input_loc_lat)
O.time = time_loop;
O.slp = output_prtmsl(:,:,nn);
O.wndspd_10m = output_wind_speed(:,:,nn);
O.wnddir = output_wind_dir(:,:,nn);
O.debug.slp = squeeze(output_prtmsl_four(:,:,nn,:));
O.debug.wndspd_10m = squeeze(output_wind_speed_four(:,:,nn,:));
O.debug.wnddir = squeeze(output_wind_dir_four(:,:,nn,:));
out_fname = sprintf('%s/HRDPS_%s.mat',fol_output,input_loc_name{nn});
save(out_fname,'-struct','O')
end

return










%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Quick plots  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf
subplot(311)
pcolor(grid_LON,grid_LAT,prtmsl)
grid on
shading flat
xlabel('Lon')
ylabel('Lat')
hold on
c = colorbar;
ylabel(c,'Pressure [Pa]')
phan=plot(input_loc_lon,input_loc_lat,'or');
legend(phan,'output points')
mycaxis = caxis;

title(datestr(time))
% Load ETOPO2 (2min global bathy) and plot coastline
Bathy = ncdataset('etopo1_bedrock_PugetSound.nc');
bathy_lat = Bathy.data('lat');
bathy_lon = Bathy.data('lon');
bathy_Band1 = Bathy.data('Band1');
contour(bathy_lon,bathy_lat,bathy_Band1,[0 0],'k')
caxis(mycaxis)

subplot(312)
pcolor(grid_LON,grid_LAT,wind_speed)
grid on
shading flat
xlabel('Lon')
ylabel('Lat')
hold on
c = colorbar;
ylabel(c,'Wind Speed [m/s]')
phan=plot(input_loc_lon,input_loc_lat,'or');
legend(phan,'output points')
mycaxis = caxis;
q_skip = 20;
quiver(grid_LON(1:q_skip:end,1:q_skip:end),grid_LAT(1:q_skip:end,1:q_skip:end),wind_u(1:q_skip:end,1:q_skip:end),wind_v(1:q_skip:end,1:q_skip:end),...
    'Color',[.5 .5 .5])

% Load ETOPO2 (2min global bathy) and plot coastline
Bathy = ncdataset('etopo1_bedrock_PugetSound.nc');
bathy_lat = Bathy.data('lat');
bathy_lon = Bathy.data('lon');
bathy_Band1 = Bathy.data('Band1');
contour(bathy_lon,bathy_lat,bathy_Band1,[0 0],'k')
caxis(mycaxis)


subplot(313)
pcolor(grid_LON,grid_LAT,wind_dir)
grid on
shading flat
xlabel('Lon')
ylabel('Lat')
hold on
c = colorbar;
ylabel(c,'Direction [deg]')
phan=plot(input_loc_lon,input_loc_lat,'or');
legend(phan,'output points')
mycaxis = caxis;
q_skip = 20;
quiver(grid_LON(1:q_skip:end,1:q_skip:end),grid_LAT(1:q_skip:end,1:q_skip:end),wind_u(1:q_skip:end,1:q_skip:end),wind_v(1:q_skip:end,1:q_skip:end),...
    'Color',[.5 .5 .5])

% Load ETOPO2 (2min global bathy) and plot coastline
Bathy = ncdataset('etopo1_bedrock_PugetSound.nc');
bathy_lat = Bathy.data('lat');
bathy_lon = Bathy.data('lon');
bathy_Band1 = Bathy.data('Band1');
contour(bathy_lon,bathy_lat,bathy_Band1,[0 0],'k')
caxis(mycaxis)
colormap(parula)

%%
% hFig = gcf;
% hFig.PaperUnits = 'inches';
% hFig.PaperSize = [8.5 13];
% hFig.PaperPosition = [0 0 8.5 13];
% print(hFig,'-dpng','-r200','Pressure_Wind_Speed')

%%
clf
pcolor(grid_LON,grid_LAT,wind_speed)
grid on
%shading flat
xlabel('Lon')
ylabel('Lat')
hold on
c = colorbar;
ylabel(c,'Wind Speed [m/s]')
phan=plot(input_loc_lon([1 16:18]),input_loc_lat([1 16:18]),'or','MarkerFaceColor','r');
legend(phan,'output points')
mycaxis = caxis;
q_skip = 5;
quiver(grid_LON(1:q_skip:end,1:q_skip:end),grid_LAT(1:q_skip:end,1:q_skip:end),wind_u(1:q_skip:end,1:q_skip:end),wind_v(1:q_skip:end,1:q_skip:end),...
    'Color',[.5 .5 .5])
xlim([-125 -121])
ylim([47 49])

% Load ETOPO2 (2min global bathy) and plot coastline
Bathy = ncdataset('etopo1_bedrock_PugetSound.nc');
bathy_lat = Bathy.data('lat');
bathy_lon = Bathy.data('lon');
bathy_Band1 = Bathy.data('Band1');
contour(bathy_lon,bathy_lat,bathy_Band1,[0 0],'k')
caxis(mycaxis)

%%
clf
pcolor(grid_LON,grid_LAT,wind_speed)
grid on
shading flat
xlabel('Lon')
ylabel('Lat')
hold on
c = colorbar;
ylabel(c,'Wind Speed [m/s]')
phan=plot(input_loc_lon,input_loc_lat,'or');
legend(phan,'output points')
mycaxis = caxis;
q_skip = 20;
quiver(grid_LON(1:q_skip:end,1:q_skip:end),grid_LAT(1:q_skip:end,1:q_skip:end),wind_u(1:q_skip:end,1:q_skip:end),wind_v(1:q_skip:end,1:q_skip:end),...
    'Color',[.5 .5 .5])


% Load ETOPO2 (2min global bathy) and plot coastline
Bathy = ncdataset('etopo1_bedrock_PugetSound.nc');
bathy_lat = Bathy.data('lat');
bathy_lon = Bathy.data('lon');
bathy_Band1 = Bathy.data('Band1');
contour(bathy_lon,bathy_lat,bathy_Band1,[0 0],'k')
caxis(mycaxis)

%%

hFig = gcf;
hFig.PaperUnits = 'inches';
hFig.PaperSize = [12 8.5];
hFig.PaperPosition = [0 0 12 8.5];
print(hFig,'-dpng','-r300','WindSpeed_Forecast_2')





%%% UNNEEDED CODE %%%%
% Create grid to interpolate on (keep approx 2.5km resolution)
% lat = 44.7:(2.5/110):56.8;
% lon = -135.2:(2.5/80):-108.6; %Note, at 56 lat, spacing is 62km, at 44 lat, 80km
