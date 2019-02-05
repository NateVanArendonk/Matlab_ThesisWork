function [ C ] = spatialcorr2( varargin)
%%Gives the Spatial Correlation between two variables (having unit time)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%spatialcorr2 is a part of TimeSaving ToolBox
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
%           (Not recommended) If you have hourly or half hourly or 6 hourly
%           data, then it is recommended to include time steps in your data
%           sets and use spatialcorr3.
%           spatialcorr2 uses interpolation of data sets on its finer grids
%           and then find the correlation of the matrix formed by
%           interpolating the data into its subgrids (Between two grid point of latitude and longitude, it divided
%           the grides into its subgrids.). Interpolating the
%           data swallow the accuracy of data sets. That's why, it''s not
%           recommended.   (For details, see the documentation)
%
%
% Syntax:
%            C=spatialcorr2(lon,lat, data1,data2 );
%
%Inputs:
%           First and second input should be longitudes and latitudes of
%           the data set. Third and fourth input should be the data between
%           which you want to find spatial correlation.
%           By default, fifth argument is set to 7. This means that
%           function divided the corresponding grids into 7 subgrids and
%           then it interpolates.
%
% Example:
%           C=spatialcorr2(lon,lat, data1,data2 );
%
%
% Please send your suggestions to the email id: ankurk017@gmail.com or
%                                               416AS2025@nitrkl.ac.in
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%spatialcorr2 is a part of TimeSaving ToolBox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin>5
    error('Sorry! Too many input arguments.')
end
dim=7;
if nargin==5
    dim=varargin{5};
end
lon=varargin{1};
lat=varargin{2};
data1=varargin{3};
data2=varargin{4};

lon1=linspace(min(lon),max(lon),(length(lon)*dim)-(dim-1));
lat1=linspace(min(lat),max(lat),(length(lat)*dim)-(dim-1));
[x,y]=meshgrid(lon1,lat1);
C1=interp2(lon,lat,data1',x',y');
C2=interp2(lon,lat,data2',x',y');
id1=(1:dim:length(lon1));
id2=(1:dim:length(lat1));
for i=1:length(id1)-1
    for j=1:length(id2)-1
        C(i,j)=nancorr2(C1(id1(i):id1(i+1)-1,id2(j):id2(j+1)-1),...
            C2(id1(i):id1(i+1)-1,id2(j):id2(j+1)-1));
    end
end
C(1,end+1)=nan;
C(end+1,:)=nan;
end
function A=nancorr2(a,b)
id1=union(find(isnan(a)),find(isnan(b)));
a(id1)=[];
b(id1)=[];
A=corr2(a,b);
end
