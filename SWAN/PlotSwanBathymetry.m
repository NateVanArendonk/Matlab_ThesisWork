
% Set file locations/names
mod_name = 'jdf';
switch mod_name
    case 'jdf'
        % Outfile name
        wind_file = 'jdf_hrdps.wnd';
        % File locations
        dep_nm = '../PS_RegionalModel/INP_sph2/jdaf_new.BOT';
        grd_nm = '../PS_RegionalModel/INP_sph2/jdf_sph_swn.grd';
        
        % From .swn file
        mxc = 659;
        myc = 1195;
        names = {'JDF','CA_JDF','Whid_Van','San_Juan_Il','Orcas_Obstruc_Blakely',...
            'Lopez_Decatur','Shaw','PennCove','Texada_Casqueti','Stuart','Spieden',...
            'Waldron','Sucia','Cyprus','Sinclair','Guemes','Vendovi','Eliza'};
        inds = [1,17,23,75,110,105,109,25,76,78,90,113,144,152,164,159,170,176];
        
    case 'ps'
        % Outfile name
        wind_file = 'ps_hrdps.wnd';
        % File locations
        dep_nm = '../PS_RegionalModel/INP_sph/PugetSound2.BOT';
        grd_nm = '../PS_RegionalModel//INP_sph/pug8_sph_swn.grd';
        
        names = {'South_PS','Whidbey','NorSeattle','Vashon','Blake','Fox',...
            'OlympiaAreaIslands','Anderson'};
        inds = [1,4,5,42,41,49,54,57];
        
        % From .swn file
        mxc = 1555;
        myc = 524;
end

% Load grid
G = swan_io_grd('read',grd_nm,mxc,myc,99,99);

% Load bottom
S.mxinp = mxc;
S.myinp = myc;
S.idla = 4;
S.nhedf = 0;
S.fname1 = dep_nm;
S.quantity = 'depth';
bot = swan_io_bot('read',grd_nm,S);

% Plot 
pcolor(G.X,G.Y,bot)
shading flat
return
%% Make Contour 
% Establish Contour Value
depth = -10;
k = [depth,depth];

% Make contour and grab the 
[C,~] = contour(G.X,G.Y,bot,k);

breaks = find(C(1,:) == depth);

% Plot Indiviual Contour Shapes
clf
for m = 1:length(breaks)
    if m ~= length(breaks) % not the last one
        lon = C(1,(breaks(m)+1):(breaks(m+1)-1));
        lat = C(2,(breaks(m)+1):(breaks(m+1)-1));
        if length(lon) > 50
            plot(lon,lat)
            disp(m)
        end
            %pause
    else
        lon = C(1,(breaks(m)+1):end);
        lat = C(2,(breaks(m)+1):end);
    end
end
% fname = sprintf('%s_%dmContour',mod_name,depth*-1);
% save(fname,'lon','lat')

%% 
for ii = 1:length(names)
    cur_name = names{ii};
    m = inds(ii); 
    lon = C(1,(breaks(m)+1):(breaks(m+1)-1));
    lat = C(2,(breaks(m)+1):(breaks(m+1)-1));
    fname = sprintf('%s_%dContour',cur_name,depth*-1);
    save(fname,'lon','lat')
    
end
%% Plot bathymetry along transect 
K.xg = zeros(length(lon));
K.yg = zeros(length(lat));
for ii = 1:length(lon)
end