dname1='y:\Projects\pscosmos\runs\psfm_4\';
fname1='salishc_b_net.nc';

info=ncinfo([dname1,fname1]);
vars=arrayfun(@(x)(x.Name),info.Variables,'un',0)';
atts=arrayfun(@(x)(x.Name),info.Variables(4).Attributes,'un',0)';
avals=arrayfun(@(x)(x.Value),info.Variables(4).Attributes,'un',0)';

mstruct=gcm(axesm('tranmerc'));
mstruct.falsenorthing=0;
mstruct.falseeasting=500000;
mstruct.geoid=[6378137 0.081819191042815];
mstruct.mapparallels=0;
mstruct.nparallels=1;
mstruct.origin=[0 -123 0];
mstruct.scalefactor=0.9996;


dname='y:\Projects\grossman\big_sac_1\';
fname='nooksack_bhb_new_trench_river_net.nc';

x=ncread([dname,fname],'NetNode_x');
y=ncread([dname,fname],'NetNode_y');
[lat,lon]=projinv(mstruct,x,y);

ncid=netcdf.open([dname,fname],'write');
varid=netcdf.inqVarID(ncid,...
    'projected_coordinate_system');
netcdf.reDef(ncid);
netcdf.renameVar(ncid,varid,'wgs84');
netcdf.close(ncid);

 cellfun(@(x,y)(ncwriteatt([dname,fname],'wgs84',...
     x,y)),atts,avals)
ncwrite([dname,fname],'NetNode_x',lon);
ncwrite([dname,fname],'NetNode_y',lat);
