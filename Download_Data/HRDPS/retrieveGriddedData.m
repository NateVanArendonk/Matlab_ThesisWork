function [ var_out ] = retrieveGriddedData( fname, fol_name, fol_location )
%retrieveGriddedData( varType, fname, fol_location )
%   By default takes second variable in gridded data file

% Set geo file
try
    nco=ncgeodataset([fol_location fol_name fname]);
    % Extract second variable, convert to double, and squeeze
    param=nco.variables(2);
    var_out=nco{param}(1,:,:);
    var_out=squeeze(double(var_out(1,:,:)));
catch
    var_out = NaN;
end


end
