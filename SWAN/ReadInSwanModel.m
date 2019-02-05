% Load in the depths for the jdf and ps model 
% ------------------------- JDF Model -------------------------------------
dep_nm = 'C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\SWAN\PS_RegionalModel\INP_sph2\jdaf_new.BOT';
grd_nm = 'C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\SWAN\PS_RegionalModel\INP_sph2\jdf_sph_swn.grd';

% From .swn file
mx = 659;
my = 1195;
% Load grid
J = swan_io_grd('read',grd_nm,mx,my,99,99);

% Load bottom
S.mxinp = mx;
S.myinp = my;
S.idla = 4;
S.nhedf = 0;
S.fname1 = dep_nm;
S.quantity = 'depth';
J.Z = swan_io_bot('read',grd_nm,S);
% -------------------------- PS Model -------------------------------------

dep_nm = 'C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\SWAN\PS_RegionalModel\INP_sph\PugetSound2.BOT';
grd_nm = 'C:\Users\ahooshmand\Desktop\PS_COSMOS\Thesis_Modeling\SWAN\PS_RegionalModel\INP_sph\pug8_sph_swn.grd';

% From .swn file
mx = 1555;
my = 524;
% Load grid
P = swan_io_grd('read',grd_nm,mx,my,99,99);

% Load bottom
S.mxinp = mx;
S.myinp = my;
S.idla = 4;
S.nhedf = 0;
S.fname1 = dep_nm;
S.quantity = 'depth';
P.Z = swan_io_bot('read',grd_nm,S);