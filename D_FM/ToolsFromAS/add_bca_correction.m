dname='y:\Projects\pscosmos\runs\psfm_m8\';
fnames=dir([dname,'*grid17_0*.cmp']);

for i=1:length(fnames)
    bca=delft3d_io_bca('read',[dname,fnames(i).name]);
    cidx=find(strcmpi('K1',bca.DATA.names));
%     bca.DATA.amp(cidx)=bca.DATA.amp(cidx)+0.0369;
     bca.DATA.phi(cidx)=bca.DATA.phi(cidx)+4.28;
    delft3d_io_bca('write',[dname,fnames(i).name],bca);
end