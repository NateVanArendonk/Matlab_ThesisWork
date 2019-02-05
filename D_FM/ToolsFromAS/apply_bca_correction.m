dname='y:\Projects\pscosmos\runs\psfm_4\';
fnames=dir([dname,'*.cmp']);

for i=1:length(fnames)
    bca=delft3d_io_bca('read',[dname,fnames(i).name]);
    bca.DATA.label='A0         0               0            ';
    sidx=find(strcmpi('K1',bca.DATA.names));
    bca.DATA.amp(sidx)=bca.DATA.amp(sidx)*1.054;
    bca.DATA.phi(sidx)=bca.DATA.phi(sidx)+6.024;
    delft3d_io_bca('write',[dname,fnames(i).name],bca);
end
