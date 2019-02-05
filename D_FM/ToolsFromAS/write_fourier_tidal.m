tstart_m=44640;
tend_m=129600;

tstart_s=tstart_m*60;
tend_s=tend_m*60;


trange=tend_s-tstart_s;
freq=[12.42;...
    12;...
    12.66;...
    23.93;...
    25.82];
fs=freq.*3600;
    


ftime=cellfun(@(x)(tstart_s+rem(trange,x):x:tend_s),...
    num2cell(fs),'un',0);
%make sure start and stop are multiples of DTUSER
dtuser=60;
fstart=cellfun(@(x)(x(1)-rem(x(1),dtuser)),ftime);
fend=cellfun(@(x)(x(end)-rem(x(end),dtuser)),ftime);

ncyc=cellfun(@(x)(numel(x)-1),ftime);

fname='y:\Projects\pscosmos\runs\psfm_m5r\fourier.fou';
fid=fopen(fname,'wt');
%write out mean first
fprintf(fid,'wl  %0.1f  %0.1f  %d  %0.3f %0.1f mean\n',...
    tstart_s,tend_s,0,1,0);
for i=1:length(fs)
    fprintf(fid,'wl  %0.1f  %0.1f  %d  %0.3f %0.1f\n',...
        fstart(i),fend(i),ncyc(i),1,0);
end
fclose(fid);
