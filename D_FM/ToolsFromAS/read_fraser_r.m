dname='f:\salishsea\data\discharge\';
fname='fraser_r.txt';
fid=fopen([dname,fname],'r');
hdr=textscan(fgetl(fid),'%s');

data=textscan(fid,'%s%f%f%f','delimiter','\t');
dstr=cell2struct(data',hdr{1});
fclose(fid);

ud=unique(dstr.DD);
mq=zeros(length(ud));
for i=1:length(ud);
    tidx=find(dstr.DD==ud(i));
    mq(i)=median(dstr.Value(tidx));
end

mins=(ud-1).*1440;

fd=fopen([dname,'fraser.tim'],'wt');
for i=1:length(mq)
    fprintf(fd,'%0.2f\t%0.2f\t%0.2f\n',...
        mins(i),mq(i),0);
end
fclose(fd);
