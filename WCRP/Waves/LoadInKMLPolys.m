% Load in polygons of WCRP locations 

kml_fol = 'C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\KML\WCRP_KML\';
kml_name = 'wcrp_locations.kml';

K = kml2struct([kml_fol kml_name]);
C = kml2struct([kml_fol 'Circles2Cut.kml']);

%% Get rid of coastal polys, only want PS
clf
for ii = 1:length(K)
    if ~inpolygon(K(ii).Lon,K(ii).Lat,C.Lon,C.Lat)
        pgon = polyshape([K(ii).Lon],[K(ii).Lat]);
        plot(pgon);
        hold on
    end
end