dname='f:\salishsea\andrew\psfm_1\';
fnames=dir([dname,'*.cmp']);

for i=1:length(fnames);
    bca=delft3d_io_bca('read',[dname,fnames(i).name]);
    bca.DATA.label='A0         1.05             0              ';
    delft3d_io_bca('write',[dname,fnames(i).name],bca);
end