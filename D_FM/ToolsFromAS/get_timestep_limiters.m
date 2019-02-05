dname='y:\Projects\pscosmos\runs\psfm_m7\DFM_OUTPUT_ps_2d\';
fname='ps_2d_map.nc';

info=ncinfo([dname,fname]);
vars=arrayfun(@(x)(x.Name),info.Variables,'un',0)';

vidx=find(strcmpi('numlimdt',vars));

nelem=info.Variables(vidx).Dimensions(1).Length;
ntime=info.Variables(vidx).Dimensions(2).Length;
limiter=ncread([dname,fname],'numlimdt',[1 ntime],...
    [nelem 1]);

lat=ncread([dname,fname],'FlowElem_ycc');
lon=ncread([dname,fname],'FlowElem_xcc');

bvals=find(limiter>0);
bmax=max(limiter);

xv=lon(bvals);
yv=lat(bvals);
vv=limiter(bvals)./sum(limiter);

dlmwrite('y:\Projects\pscosmos\runs\psfm_m7\time_step_limiters.xyz',...
    [xv yv vv],'delimiter','\t','precision','%0.6f')