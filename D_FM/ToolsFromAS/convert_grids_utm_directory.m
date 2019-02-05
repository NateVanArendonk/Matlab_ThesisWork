dname='z:\projects\pscosmos\model\andrew\';
fnames=dir([dname,'*.grd']);
grd=arrayfun(@(x)(wlgrid('read',[dname,x.name])),fnames);

mstruct=gcm(axesm('tranmerc'));
mstruct.falsenorthing=0;
mstruct.falseeasting=500000;
mstruct.geoid=[6378137 0.081819191042815];
mstruct.mapparallels=0;
mstruct.nparallels=1;
mstruct.origin=[0 -123 0];
mstruct.scalefactor=0.9996;
[x,y]=arrayfun(@(x)(projfwd(mstruct,x.Y,x.X)),grd,'un',0);
[grd(:).X]=deal(x{:});
[grd(:).Y]=deal(y{:});

[grd(:).CoordinateSystem]=deal('Cartesian');
[~,fp]=arrayfun(@(x)(fileparts(x.FileName)),grd,'un',0)
fn=cellfun(@(x)([x,'_utm.grd']),fp,'un',0)
[grd(:).FileName]=deal(fn{:});

arrayfun(@(x)(wlgrid('write',[dname,x.FileName],x)),grd)