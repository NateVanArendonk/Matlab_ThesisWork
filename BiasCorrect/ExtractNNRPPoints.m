% This code will extract individual NNRP Model Points
% It is not automated, it will require user input to pick each point of
% interest and save that variable to the matlab workspace 
clearvars
N = load('E:\Abbas\PS_COSMOS\Thesis_Modeling\Quantile_Correction\WA_NNRP.mat');

% Grab land and water points
land_lat = N.lat_m(N.land_mask);
land_lon = N.lon_m(N.land_mask);
water_lat = N.lat_m(~N.land_mask);
water_lon = N.lon_m(~N.land_mask);

% Load Coastline
load('E:\Abbas\PS_COSMOS\Thesis_Modeling\Quantile_Correction\WA_coast.mat');
% Find all shapes above threshold
wa_lat = [];
wa_lon = [];
thresh = 500 ;
for j = 1:length(wa_coast)
    temp_x = wa_coast(j).X;
    temp_y = wa_coast(j).Y;
    if length(temp_x) >= thresh %&& j ~= 3348 % 3348 is oregon
        for m = 1:length(temp_x)
            wa_lat(end+1) = temp_y(m);
            wa_lon(end+1) = temp_x(m);
        end
    end
end

% Plot model points
clf
plot(wa_lon,wa_lat,'k')
hold on
l = plot(land_lon,land_lat,'.r','MarkerSize',14);
w = plot(water_lon,water_lat,'b.','MarkerSize',14);
xlim([-123.2 -122])
ylim([47.1 48])

%% Pick NNRP points and extract them to matlab workspace and save 


% User picks point and exports data to workspace then changes name and
% presses play
name = 'swan_69_NNRP_Point';
plon = p.Position(1);
plat = p.Position(2);
save(name,'plon','plat');

fname = sprintf('%s.mat',name);
% movefile(fname,'NNRP_Points')


%% Make Plot of extracted points
clf
plot(wa_lon,wa_lat,'k')
hold on
axis equal
xlim([-122.8 -122.2])
ylim([47 47.7])
xlabel('Degrees Longitude')
ylabel('Degrees Latitude')
set(gca,'FontSize',14)

np = dir('NNRP_Points\*.mat');
for ii = 1:length(np)
    Nstring = strfind(np(ii).name,'N');
    Pstring = strfind(np(ii).name,'.');
    num = np(ii).name(Nstring+1:Pstring-1);
    np(ii).num = str2num(num);
end
np = nestedSortStruct(np,{'num'});
b = parula(length(np));
Z = struct;
for f = 1:length(np)
%     f_open = strcat('NNRP_Points\',np(f).name);
    n = load(strcat('NNRP_Points\',np(f).name));
    [~,Z(f).name,~] = fileparts(strcat('NNRP_Points\',np(f).name));
    ll(f) = plot(n.plon,n.plat,'o','Color',b(f,:),'MarkerFaceColor',b(f,:),'MarkerSize',12);
end
% l = plot(land_lon,land_lat,'.r','MarkerSize',14);
% w = plot(water_lon,water_lat,'k.','MarkerSize',14);
% ll(19) = l;
% ll(20) = w;
% Z(19).name = 'NNRP Land';
% Z(20).name = 'NNRP Water';

% lgd = legend(ll,arrayfun(@(x)(x.name),Z,'un',0),...
%     'interpreter','none',...
%     'location','northwest');
% lgd.FontSize = 12;
% printFig(gcf,'TacomaExtractedPointsNNRP',[15 13],'png',300)

%% Save NNRP Timeseries for each point in NNRP Points folder
np = dir('NNRP_Points\*.mat'); % List of NNRP Points that I want to save
for f = 1:length(np)
    fol_name = 'temp\';
    if ~exist([fol_name np(f).name])
        fprintf('Saving %s\n',np(f).name)
        f_open = strcat('NNRP_Points\',np(f).name);
        [~,name,~] = fileparts(f_open);
        n = load(f_open);
        % Find position on NNRP Grid
        if isfield(n,'plat')
            [r,c] = find(N.lat_m == n.plat & N.lon_m == n.plon);            
        else
            [r,c] = find(N.lat_m == n.Position(2) & N.lon_m == n.Position(1));
        end
        % Grab data at exact location
        elev = N.elev_m(r,c);
        slp = squeeze(N.slp(r,c,:));
        wndspd = squeeze(N.wndspd(r,c,:));
        wnddir = squeeze(N.wnddir(r,c,:));
        surfp = squeeze(N.surfp(r,c,:));
        time = N.time;
        if isfield(n,'plat')
            lat = n.plat;
            lon = n.plon;
        else
            lat = n.Position(2);
            lon = n.Position(1);
        end
        land = land_mask(r,c);
        % Save the data
        outname = sprintf('%s.mat',name);
        save(outname,'elev','slp','wndspd','wnddir','surfp','time','lat','lon','land')
        movefile(outname,'temp')
    end
end
