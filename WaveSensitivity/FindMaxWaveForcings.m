%Creates a hindcast for a specific transect and calculates the Forcings for
%Max Hsig for storm characterization
%   This script will create a hindcast along any specific transect created
%   in google earth.  The user determines what type of data they want to
%   use as a forcing.  The data can be obs, quantile corrected NNRP, or
%   NNRP.  This function should be called in a script and then the output
%   should be used to plot.  



clearvars
addpath C:\Functions_Matlab
addpath C:\Functions_Matlab\cbrewer
addpath C:\Functions_Matlab\time
addpath C:\Functions_Matlab\ocean
addpath C:\Functions_Matlab\mapping\kml

% Load Data Below


slr_val = 2;
m2ft = .3048;
slr = slr_val*m2ft; %[m]

% Folder Locations 
nnrp_fol = 'C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\Quantile_Correction\NNRP_PointData\';
mask_fol = 'C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\KML\Basin_Masks\';
kml_fol = 'C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\KML\PugetSoundAreas\';
lut_fol = 'C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\Hindcast\LUT_Extracts\';
% kml_fol = 'C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\KML\PugetSoundShorline\';
% 
% % Make key-value pairs for nnrp_points and polygons 
% keys = {'HoodCanal.kml','NorthCentral.kml','Olympia.kml','OlympiaNorth.kml',...
%     'OlympiaWest.kml','PortTownsend_Hood.kml','Skagit.kml','SouthCentral.kml',...
%     'Tacoma.kml','Whidbey_JDF.kml'};
% values = {'Nor_Hood_Canal_NNRP_Point.mat','Kingston_NNRP_Point.mat','South_Tacoma_NNRP_Point.mat',...
%     'Nor_Olympia_NNRP_Point.mat','Middle_Olympia_NNRP_Point.mat','Townsend_Whidbey_NNRP_Point.mat',...
%     'Everett_NNRP_Point.mat','Blake_Island_NNRP_Point.mat','Tacoma_NNRP_Point.mat','JDF_Whidbey_NNRP_Point.mat'};
% M = containers.Map(keys,values); % Creates a dictionary of key value pairs 

% --- Load LUT 
lut2load = [lut_fol 'seatac_LUT_offshore_extract.mat'];
L = load(lut2load);
[L.x_utm,L.y_utm,L.utmzone] = deg2utm(L.y,L.x);
% --- Load in Tides
load('C:\Users\ahooshmand\Desktop\Data_Forcings\COOPS_tides\seattle\seattle_hrV_NAVD88.mat');
T.twl = tides.WL_VALUE;
T.time = tides.time;

% --- Load KML
% K = kml2struct([kml_fol 'wa_10m_FULL.kml']);
K = kml2struct([kml_fol 'Seattle_Tacoma.kml']);

% --- Load in single NNRP Point for the time vector 
nn2load = [nnrp_fol 'blakeisland_corwith_westpoint.mat'];
N = load(nn2load);
% Fix Time
N.time  = N.time';
time = N.time;
% clear N
disp('Done Loading Data')
return
%% Create Hindcast from LUT

% Interpolate twl onto time vec
inds = find(unique(T.time));
T.time = T.time(inds);
T.twl = T.twl(inds);
inds = find(diff(T.time) == 0);
T.time(inds) = [];
T.twl(inds) = [];
T.twl_i = interp1(T.time,T.twl,time);

% Dimensions of LUT
twl = T.twl_i; twl = double(twl)';
% setup for interp
[X,Y,Z] = meshgrid(L.tide,L.speed,L.dir);

% Loop over extraction points
tic
% Initialize
% ---1. Yearly Max
S.date = NaN(length(L.x),1);
S.spd = S.date;
S.dir = S.date;
S.wl = S.date;
S.hs = S.date;
S.tp = S.date;
S.ind = S.date;

    speed = squeeze(N.wndspd);
    wnddir = squeeze(N.wnddir);
for nn = 1:length(L.x)
    if mod(nn,100) == 0
        fprintf('%d out of %d\n',nn,length(L.x))
    end
% %     fprintf('%d out of %d\n',nn,length(L.x))
%     % Find which forcing to use based on point being in mssk polygon
%     masks = dir('C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\KML\Basin_Masks\*.kml'); % Grab all of the masks 
%     for f = 1:length(masks) % loop through the masks 
%         f_open = strcat(mask_fol,masks(f).name);
%         [~,name,ext] = fileparts(f_open);
%         m = kml2struct(f_open);
%         % See if point is inside polygon 
%         in = inpolygon(L.x(nn),L.y(nn),m.Lon,m.Lat);
%         if in 
%             break
%         end
%         
%     end
    
    % Grab NNRP Forcing 
%     NN_point = M([name ext]);
%     N = load([nnrp_fol NN_point]);
%     speed = squeeze(N.wndspd);
%     wnddir = squeeze(N.wnddir);
    
    % Grab LUT data and determine Wave stats based on NNRP
    %Hsig
    temp = permute(squeeze(L.hs(nn,:,:,:)),[2 1 3]); 
    hs_ts = interp3(X,Y,Z,temp,twl,speed,wnddir,'spline')'; % Hsig timeseries
%     hs_ts_slr = interp3(X,Y,Z,temp,twl+slr,speed,wnddir,'spline')';

    % Period
    temp = permute(squeeze(L.tp(nn,:,:,:)),[2 1 3]);
    tp_ts = interp3(X,Y,Z,temp,twl,speed,wnddir,'spline')';
%     tp_ts_slr = interp3(X,Y,Z,temp,twl+slr,speed,wnddir,'spline')';
    
    % Calculate Metrics about forcings
    [Mx,I] = max(hs_ts);
    S.date(nn) = time(I);
    S.spd(nn) = speed(I);
    S.dir(nn) = wnddir(I);
    S.wl(nn) = twl(I);
    S.hs(nn) = Mx;
    S.tp(nn) = tp_ts(I);
    S.ind(nn) = I;
end
toc
disp('Done Creating Hindcast')
clear N T X Y Z twl tides temp speed dir
S.x = L.x_utm;
S.y = L.y_utm;
S.lat = L.y;
S.lon = L.x;
S.time = time;
%save
save('MaxWaveForcings_SeattleTacoma','-struct','S')
disp('Transect Saved')
toc
return
clearvars
load('MaxWaveForcings_SeattleTacoma.mat')
return
%% Split storms up and grab indices 
[B,I] = sort(date);
s = SplitVec(B);
inds2grab = 0;
% Split up storms 
for ii = 1:length(s)
    S(ii).date = s{ii,1}(1); % Date of event 
    [r,c] = size(s{ii,1});
    if ii == 1   
        S(ii).cells = {I(1:r)}; % Location along transect 
        S(ii).spd = spd(I(1:r)); % Speed of forcing
        S(ii).dir = dir(I(1:r)); % Direction of forcing
        S(ii).wl = wl(I(1)); % WL at time of event
        S(ii).storm_inds = length(S(ii).spd);
        inds2grab = inds2grab + r;
        
    else
        S(ii).cells = {I(inds2grab+1:inds2grab+r)};
        S(ii).spd = spd(I(inds2grab+1:inds2grab+r));
        S(ii).dir = dir(I(inds2grab+1:inds2grab+r));
        S(ii).wl = wl(I(inds2grab+1));
        S(ii).storm_inds = length(S(ii).spd);
        inds2grab = inds2grab+r;
    end
end

% Sort storms from most to least number of cells
storm_length=arrayfun(@(x)(length(x.spd)),S); 
[s_sort,sidx]=sort(storm_length,'descend');
S2 = S(sidx);
S = S2;
clear S2

% Change speed and direction to only 1 value - only do this for singel
% domain analysis!
for ii = 1:length(S)
    S(ii).spd = S(ii).spd(1);
    S(ii).dir = S(ii).dir(1);
    S(ii).percent = 100*(S(ii).storm_inds/length(B));
end
%% Plot storm results along shore
clf
IM = load('C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\GeoTiffs/seattle_tacoma3.mat');
imagesc(IM.xm,IM.ym,IM.im);
set(gca,'ydir','normal');
hold on

b = parula(length(S));
for ii = 1:length(S)
    S(ii).name = sprintf('Storm %d',ii);
    S(ii).fromDirection = wrap2360(S(ii).dir +180);
    plot(x(S(ii).cells{1,1}),y(S(ii).cells{1,1}),'o','MarkerSize',9,'MarkerFaceColor',b(ii,:),'MarkerEdgeColor',b(ii,:))
    hold on
end
lgd = legend(arrayfun(@(x)(num2str(x.name,3)),S,'un',0),...
    'interpreter','none',...
    'location','northwest'); % Legend Function written by Andrew Stevens
% title(lgd,'Direction of Max Hsig')
axis equal
xlim([5.30*10^5 5.55*10^5]);
ylim([5.232*10^6 5.28*10^6]);
xlabel('Easting [m]')
ylabel('Northing [m]')
printFig(gcf, 'SeattleTacoma_WindDirectionMaxWave',[8.5 11],'png',300)
%% Write CSV 
fid = fopen('SeattleTacoma_StormWaves.csv','w');

% Write data to new file
fprintf(fid,'Datenum, Date, Speed, Direction, Water Level, Number of Indices, Percent of Shoreline Impacted\n');
for ii = 1:length(S);
    fprintf(fid,'%.2f, %s, %.2f, %.2f, %.2f, %d, %.2f\n',S(ii).date, datestr(S(ii).date), S(ii).spd, S(ii).dir, S(ii).wl, S(ii).storm_inds, S(ii).percent);
end
    
fclose(fid)