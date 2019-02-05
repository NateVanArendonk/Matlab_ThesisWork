function [gx,gy] = getXYgridParams(xygrid)
%Given the OET created .grd file, extract the size and subtract one, needed
%for input to Xbeach
%   Outputs: # in x direction(gx) and in y direction (gy)
%   Input: Full path name of .grd file 

fid = fopen(xygrid); % Open .grd file
format = '%s %s'; % Format for reading in 
data = textscan(fid,format,'HeaderLines',1); % 
fclose(fid);

tx = str2double(data{1}{1});
ty = str2double(data{2}{1});

gx = tx - 1;
gy = ty - 1;
end

