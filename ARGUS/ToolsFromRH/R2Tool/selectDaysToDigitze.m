% find wave conditions for months of Aug through Nov, 2015
dn1 = datenum(2015,8,1,17,0,0);
dn2 = datenum(2015,11,30,17,0,0);
dns = dn1:dn2;

for i = 1: length(dns)
    [zt(i), Hs(i), fp(i)] = getDuckEnvirInfo(matlab2Epoch(dns(i)));
end

figure(1); clf
plot(dns, Hs)
datetick('x'); grid on

N = 15;
[HSort,ind] = sort(Hs);
pick = (5:round(length(dns)/N): length(dns));
HSort(pick)
dates = dns(ind(pick));
datestr(dates)

NRMax = [52 64 67 9 88 71 72 73 79 35 36 37 38 48 49 50 51 53 54 55];        % days to digitize from RMax
HsForDays = Hs(NRMax)
dnForDays = dns(NRMax)