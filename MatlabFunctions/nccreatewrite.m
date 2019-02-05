function C= nccreatewrite(varargin)
%%Creates the netcdf4 file directly.
%%Get rid of writing the long commands everytime to create the netcdf4 file
%%or to add the variable to the nc file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nccreatewrite is a part of TimeSaving ToolBox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Created by Ankur Kumar                        Wednesday; January 03, 2018
%                                              Version: 1.0
%
% National Institute of Techonology, Rourkela
% Rourkela, Odisa - 769008, India
% Department of Earth and Atmospheric Sciences
% Email: ankurk017@gmail.com
%        416AS2025@nitrkl.ac.in
%
%
% Function:
%           One can create a netcdf4 file with the predefined MATLAB
%           function using nccreate and ncwrite. This function is a
%           combination of these two predefined functions. One can use this
%           to save the time and get rid of writing the same commands for
%           storing multi variables in a nc file.
%
% Syntax:
%           nccreatewrite('test1.nc','lat',{'lat','c'},lat)
%
%Inputs:
%           First input should be the name of the nc file in which you want
%           to store the data
%           Second input should be the name of the variable in which you
%           want to store the specific variable
%           Third argument should be in braces (not structure) which must
%           contains the number of variables as that of the size of the
%           data you want to store.
%           ex: If you want to store 'lat' whose dimensions is 5*1, then
%           third argument should be {'a','b'}
%           The above is becasue of MATLAB also stores the variable
%           dimensions seprately, not in the variable list.
%           You can see this when you use ncdisp to see the listed
%           variables in nc file.
%           Fourth argument should be the data you want to write in nc
%           file.
%
% Example:
%            clc
%            clear
%            lon=(65:0.5:95)';
%            lat=(3:0.5:35)';
%            data=randi(20,65,61,365);
%            delete test1.nc
%            nccreatewrite('test1.nc','lat',{'lat','c'},lat)
%            nccreatewrite('test1.nc','lon',{'lon','c'},lon)
%            nccreatewrite('test1.nc','TC',{'lat','lon','days'},data)
%
% 
% Please send your suggestions to the email id: ankurk017@gmail.com or
%                                               416AS2025@nitrkl.ac.in
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nccreatewrite is a part of TimeSaving ToolBox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin>=5
    error('Only 4 argument needed. Type help nccreatewrite to get help of this function')
end

filename=varargin{1};
variable_name=varargin{2};
var=varargin{3};
data=varargin{4};

AA=size(data);
A=length(AA);
B=length(var);
if A~=B
    error('Dimensions of the data which you want to write is not matching with the given number of arguments for storing the dimensions of data.')
end
CC=[];
for i=1:A
    C={var{i},AA(i)};
    CC=[CC C];
end
nccreate(filename,variable_name,...
    'Dimensions',CC,...
    'format','netcdf4')
ncwrite(filename,variable_name,data);
end