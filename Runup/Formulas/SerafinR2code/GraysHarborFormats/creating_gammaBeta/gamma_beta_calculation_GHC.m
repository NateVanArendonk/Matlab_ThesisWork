cd C:\data\projects\GraysHarbor
load('wash_2002.mat')
wash_2002.x = wash_2002.x./1000;
wash_2002.y = wash_2002.y./1000;
indy = find(wash_2002.lat>=46.7 & wash_2002.lat<=47.54); %GH COUNTY



figure(1)
plot(wash_2002.x(indy),wash_2002.y(indy),'b.')
hold on

% Interpolate a line?
% First break the coast up into 1 km bins?
bins = 5170.5:1:5270.5;
hline(bins,'k')
% FIND ALL POINTS IN BINS AND CREATE A LINE BETWEEN THOSE POINTS
x = wash_2002.x(indy);
y = wash_2002.y(indy);
clear m

for ii = 1:length(bins)-1
    ind = find(y>=bins(ii)&y<=bins(ii+1));
    if isempty(ind)
        b(ii) =NaN;
        m(ii) = NaN;
        east_v(ii) = NaN;
        west_v(ii) = NaN;
    else
        xr = [ones(length(ind),1) x(ind)];
        B = regress(y(ind),xr);
        %slope = (y(ind(end))-y(ind(1)))/(x(ind(end))-x(ind(1)));
        yr = B(1)+B(2).*x(ind);
        %bb = y(ind(1))-(slope*x(ind(1)));
        %yr = bb +slope.*x(ind);
        plot(x(ind),yr,'g','linewidth',2)
        b(ii) = B(1);
        m(ii) = B(2);
        %b(ii) = bb;
        %m(ii) = slope;
        east_v(ii) = nanmean(x(ind));
        west_v(ii) = nanmean(y(ind));
    end

end

% First create evenly spaced points throughout the dataset
dd = atand(m);
newdeg = dd;
% How far from 90? deal with negative numbers first

ind = find(dd<0);
oops = 90+dd(ind);
negVals = 270-oops;
newdeg(ind) = negVals;
% How far from 90? deal with positive numbers second!
ind = find(dd>0);
oops = 90-dd(ind);
posVals = 270+oops;
newdeg(ind) = posVals;

figure(4)
plot(newdeg,bins(1:end-1),'g.-')

figure
scatter(east_v,west_v,50,newdeg,'filled')
colorbar
grid on
set(gca,'linewidth',2,'fontsize',14','fontweight','demi')
xlabel('Easting (km)')
ylabel('Northing(km)')
print('shorelineangle_bins.png','-dpng','-r300')

[lat, lon] = utm2ll(west_v*1000, east_v*1000, 10);
gamma_b.x = east_v*1000;
gamma_b.y = west_v*1000; 
gamma_b.lat = lat;
gamma_b.lon = lon;
gamma_b.val = newdeg;


gamma_beta = nan(length(y),1);
for ii = 1:length(bins)-1
    ind = find(y>=bins(ii)&y<=bins(ii+1));
    if isempty(ind)
    else
    gamma_beta(ind) = newdeg(ii);
    end

end

figure
scatter(x,y,50,gamma_beta,'filled')
colorbar
grid on
set(gca,'linewidth',2,'fontsize',14','fontweight','demi')
xlabel('Easting (km)')
ylabel('Northing(km)')
print('shorelineangle.png','-dpng','-r300')

