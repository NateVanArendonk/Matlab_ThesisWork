dname='z:\projects\pscosmos\data\bathy\';
xyzfiles=dir([dname,'*.xyz']);
ascfiles=dir([dname,'*.asc']);

xyz1=arrayfun(@(x)(load([dname,x.name])),xyzfiles,'un',0);
xyz2=arrayfun(@(x)(ascii2xyz([dname,x.name])),ascfiles,'un',0);

xyz=cat(1,xyz1,xyz2);
clear xyz1 xyz2



%get the global limits
fun=@(x)([min(x) max(x)]);
xlims=fun(xyz{1}(:,1));
ylims=fun(xyz{1}(:,2));


xdist=vdist(ylims(1),xlims(1),ylims(1),xlims(2));
xi=linspace(xlims(1),xlims(2),floor(xdist/100));
ydist=vdist(ylims(1),xlims(1),ylims(2),xlims(1));
yi=linspace(ylims(1),ylims(2),floor(ydist/100));

[X,Y]=meshgrid(xi,yi);

T=scatteredInterpolant(xyz{1}(:,1),xyz{1}(:,2),...
    xyz{1}(:,3),'nearest');
gb=T(X,Y);

rank=[4 2 1 3];
Z=bin2mat(xyz{rank==1}(:,1),xyz{rank==1}(:,2),...
    xyz{rank==1}(:,3),X,Y);
for i=2:length(rank)-1
    zt=bin2mat(xyz{rank==i}(:,1),xyz{rank==i}(:,2),...
        xyz{rank==i}(:,3),X,Y);
    zmask=isnan(Z);
    Z(zmask)=zt(zmask);
end
Z(isnan(Z))=gb(isnan(Z));

cs=contour(X,Y,Z,[1 1]);
cdata=contourdata(cs);


dout='z:\projects\pscosmos\data\bathy\';
fout='salish_sea_1m.ldb';
fid=fopen([dout,fout],'wt');
for i=1:length(cdata)
    fprintf(fid,'C%3.3d\n',cdata(i).level);
    fprintf(fid,'%d %d\n',cdata(i).numel,2);
    
    for j=1:cdata(i).numel
        fprintf(fid,'%0.6f  %0.6f\n',cdata(i).xdata(j),...
            cdata(i).ydata(j));
    end
end
fclose(fid);





