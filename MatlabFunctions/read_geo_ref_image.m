function [ xim, yim, im ] = read_geo_ref_image( image_file, ref_file )
%[ xim, yim, im ] = read_geo_ref_image( image_file, ref_file )
%   Works with tif/tfw, png/pgw, etc 

im=imread(image_file);
r=worldfileread(ref_file);
[xim,yim]=pixcenters(r,size(im));

end

