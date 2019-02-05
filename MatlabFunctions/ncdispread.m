function [ Z,aa ] = ncdispread( varargin )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ncdispread is a part of TimeSaving ToolBox
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
%           ncdispread displays and store all the variables at a time. No
%           need to load the variables seprately.
% Syntax:
%           [A a]=ncdispread('sample_nc_file.nc')
%
%Inputs:
%           It takes only one input as the name of the file. If you store
%           the output in one variable, then it only displays and loads all datas of all
%           variables in the cell foramt but you can not see the names of each variables. So,
%           it is recommened to store output in two variables.
%
%
% Example:
%
%            [A a]=ncdispread('sample_nc_file.nc')
%
%
% Please send your suggestions to the email id: ankurk017@gmail.com or
%                                               416AS2025@nitrkl.ac.in
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ncdispread is a part of TimeSaving ToolBox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin>=2
    error('Only one argument is required. Enter the name of nc file to display and load all the variables.')
end
nc_file_name=varargin{1};
ncdisp(nc_file_name)

a=dir(nc_file_name)
size=((a.bytes)/(1e6));

if size>1500
    prompt = 'File Size is more than 1.5GB. Loading all the variables at a time takes lot of time and may be your MATLAB freeze up. \n Do you still want more? Y/N [Y]: ';
    str = input(prompt,'s');
    if str == 'N'
        return
    end
end
A=ncinfo(nc_file_name);
for i=1:length(A.Variables)
    aa{i}=A.Variables(i).Name;
    Z{i}=ncread(nc_file_name,aa{i});
end
end



