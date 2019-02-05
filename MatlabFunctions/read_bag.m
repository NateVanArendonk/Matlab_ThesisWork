function bag = read_bag(varargin)
% READ_BAG - reads data from a bathymetric attributed grid file.
%   BAG = read_bag(BAGFILE) reads the contents of a bathymetric 
%   attributed grid file BAGFILE (string input). Optional inputs
%   described below facilitate reading a subset of the gridded 
%   dataset in the BAGFILE.
% 
% OPTIONAL INPUTS:
%   BAG = read_bag(...,'xlim',[minx maxx],'ylim',[miny maxy]) 
%      Optional input arguments 'xlim' and 'ylim' can be supplied to 
%      limit the spatial extents of the grid.
%
%   BAG = read_bag(...,'stride',X) - The keyword 'stride' can be 
%      supplied along with a scalar to decrease the resolution of
%      the grid output, and READ_BAG will read every X value in the grid. 
%      The value of 'stride' must be a positive integer. 
%
%   BAG = read_bag(...,'read_unc','true') - Reads uncertainty grid
%       as well as elevations.  Default = 'false'. 
%
% OUTPUT: 
%   The contents of the bag file are output as a structure array 
%   BAG with the following fields:
%       'filename'     - Bag file name
%       'projection'   - Coordinate system projection specified in 
%                        bag file metadata
%       'zone'         - Zone designator for projection
%       'datum'        - Horizontal datum 
%       'vert_datum'   - Vertical datum
%       'grid_extents' - Corner points of grid ([x1 x2 y1 y2])
%       'x_resolution' - Resolution of output grid in x direction
%       'y_resolution' - Resolution of output grid in y direction
%       'x'            - X coordinates (1 x n) of grid
%       'y'            - Y coordinates (m x 1) of grid
%       'z'            - Z values (m x n) in grid
%       'unc'          - Uncertainty (m x n) in grid values (optional)
% 
%   Currently only a portion of the metadata included in the bag file is
%   read. The .bag format is relies on a hdf5 file structure and this code
%   uses the built-in H5READ function (introduced in R2011a). More 
%   information about the bag file format can be found here:
%   https://marinemetadata.org/references/bag
% 
% SEE ALSO h5read

% Andrew Stevens @ USGS
% astevens@usgs.gov
% 8/31/2015

%parse inputs
narginchk(1,inf)

p=inputParser;
p.KeepUnmatched=1;
opts={'filename',[],      {'char'},       {};...
    'xlim',      [],      {'numeric'},    {'numel',2};...
    'ylim',      [],      {'numeric'},    {'numel',2};...
    'stride',    1,       {'numeric'},    {'integer';'positive'};...
    'read_unc',  'false', {'char'},       {}};

cellfun(@(x)(p.addOptional(x{1},x{2},...
    @(y)(validateattributes(y, x{3},x{4})))),num2cell(opts,2))
p.parse(varargin{:});
opt=p.Results;
validatestring(opt.read_unc,{'true';'false'},'read_bag','read_unc');

%make sure it exists
if ~exist(opt.filename,'file')
    error('File not found.')
end


%need to parse metadata in a xml struct
info = h5info(opt.filename,'/BAG_root/metadata');
mdata = h5read(opt.filename,'/BAG_root/metadata');

mstr=cat(2,mdata{:});
if ~isempty(info.FillValue)
    idx=strfind(mstr,info.FillValue);
    mstr(idx)=[];
end


stream = java.io.StringBufferInputStream(mstr);
factory = javaMethod('newInstance', ...
    'javax.xml.parsers.DocumentBuilderFactory');
builder = factory.newDocumentBuilder;
doc = builder.parse(stream);

%collect some of the metadata..could be expanded later
fields={'title';...
    'resolution';...
    'cornerPoints';...
    'projection';...
    'zone';...
    'ellipsoid';...
    'datum'};
vals=cell(size(fields,1),1);
%running into problems with xml parsing for many files

for i=1:length(fields)
    list=doc.getElementsByTagName(fields{i});
    len = list.getLength;
    vals{i} = cell(1, len);
    for j=1:len
        vals{i}(j) = list.item(j-1).getTextContent;
    end
end
    
    

%extract the corner points
fmt='%f,%f %f,%f';
cp=cell2mat(textscan(char(vals{3}),fmt));

%figure out which grid points to read
xres=str2double(vals{2}(2));
yres=str2double(vals{2}(1));
xall=cp(1):xres:cp(3);
yall=cp(2):yres:cp(4);

fun=@(x,y)(find(x>=y(1) & x<=y(2)));
%x-axis
if ~isempty(opt.xlim)
    xv=fun(xall,opt.xlim);
    if isempty(xv)
        error(['No data found within specified limits. \n',...
            'Grid extents are : %f %f %f %f'],...
            cp(1),cp(3),cp(2), cp(4));
    else
        xcount=floor(numel(xv)/opt.stride);
    end
else
    xcount=floor(numel(xall)/opt.stride);
    xv=1;
    
end


%y-axis
if ~isempty(opt.ylim)
    yv=fun(yall,opt.ylim);
    if isempty(yv)
        error(['No data found within specified limits. \n',...
            'Grid extents are : %f %f %f %f'],...
            cp(1),cp(3),cp(2), cp(4));
    else
        ycount=floor(numel(yv)/opt.stride);
    end
else
    ycount=floor(numel(yall)/opt.stride);
    yv=1;
end

xg=xall(xv(1):opt.stride:xv(1)+((xcount-1)*opt.stride));
yg=yall(yv(1):opt.stride:yv(1)+((ycount-1)*opt.stride));

%package the metadata
bag=struct('filename',vals{1},...
    'projection',vals{4},...
    'zone',vals{5},...
    'datum',vals{6},...
    'vert_datum',vals{7}(2),...
    'grid_extents',[min(xg) max(xg) min(yg) max(yg)],...
    'y_resolution',yres*opt.stride,...
    'x_resolution',xres*opt.stride);

%generate x and y vals and read dataset
bag.x=xg;
bag.y=yg';
bag.z=double(h5read(opt.filename,'/BAG_root/elevation',...
    [xv(1) yv(1)],[xcount ycount],[opt.stride opt.stride])');
minz = h5readatt(opt.filename,...
    '/BAG_root/elevation','Minimum Elevation Value');
maxz = h5readatt(opt.filename,...
    '/BAG_root/elevation','Maximum Elevation Value');
bag.z(bag.z<minz | bag.z>maxz)=NaN;


if strcmp(opt.read_unc,'true')
    bag.unc=double(h5read(opt.filename,'/BAG_root/uncertainty',...
        [xv(1) yv(1)],[xcount ycount],[opt.stride opt.stride])');
    
    minu = h5readatt(opt.filename,...
        '/BAG_root/uncertainty','Minimum Uncertainty Value');
    maxu = h5readatt(opt.filename,...
        '/BAG_root/uncertainty','Maximum Uncertainty Value');
    bag.unc(bag.unc<minu | bag.unc>maxu)=NaN;
end
    

