function matlab2kml(name,lonlat)
%makes a kml file for use in google earth from a matrix of lat lon values
%input:  name of kml file, one matrix containing longitude and latitude
%variables.
%Name is name of output kml file 
% Make sure lonlat variable is a #x2 matrix with Lon first then Lat
%usage:  matlab2kml('wa_coastline',latlon)

% Adapted from pwr_kml function from matlab file exchange.  

% Headers to be used
header1 = '<?xml version="1.0" encoding="UTF-8"?>';
header2 = '<kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http://www.google.com/kml/ext/2.2" xmlns:kml="http://www.opengis.net/kml/2.2" xmlns:atom="http://www.w3.org/2005/Atom">';

fid = fopen([name '.kml'], 'wt');
d=flipud(rot90(fliplr(lonlat)));
fprintf(fid, '%s \n',header1);
fprintf(fid, '%s \n',header2);
fprintf(fid, '<Document> \n');
fprintf(fid, '        <Placemark> \n');
fprintf(fid, '                <name>%s</name> \n',name);
fprintf(fid, '                <LineString> \n');
fprintf(fid, '                        <tessellate>1</tessellate> \n');
fprintf(fid, '                        <coordinates> \n');
for nn = 1:length(lonlat)
    if nn == 1
        fprintf(fid, '                                %.6f,%.6f,0 ', lonlat(nn,1),lonlat(nn,2));
    elseif nn == length(lonlat)
        fprintf(fid, '%.6f,%.6f,0 \n',lonlat(nn,1),lonlat(nn,2));
    else
        fprintf(fid, '%.6f,%.6f,0 ',lonlat(nn,1),lonlat(nn,2));
    end
end
fprintf(fid, '                        </coordinates> \n');
fprintf(fid, '                </LineString>\n');
fprintf(fid, '         </Placemark>\n');
fprintf(fid, '</Document>\n');
fprintf(fid, '</kml>');
fclose(fid)