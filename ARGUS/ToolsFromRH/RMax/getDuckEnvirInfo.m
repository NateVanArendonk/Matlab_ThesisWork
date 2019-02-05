function [ztm, Hs, fp] = getDuckEnvirInfo(et)
%   [ztm, Hs, fp] = getDuckEnvirInfo(eTime)
%
% Finds the measured vertical tide elevation, ztm, at time eTime.  
% Also find the significant wave height, Hs, and a
% peak wave frequency, fp.  All times are epoch.  
% ztm and Hs are in meters, fp is in Hz.

pnTides = '/ftp/pub/argus02b/tides/';
pnWaves = '/ftp/pub/argus02b/waves/';
tidefNames = 'FRF-ocean_waterlevel_eopNoaaTide_2015';
wavesfNames = 'FRF-ocean_waves_waverider-17m_2015';

[yr mon dy hr mn sec] = datevec(epoch2Matlab(et));
if yr~=2015
    error('Only 2015 has current data')
end

in = [pnTides tidefNames num2str(mon, '%02d') '.nc'];
outTide = parseNCFile(in);
t = outTide.var.time;
zt = outTide.var.waterLevel;
ztm = interp1(t,zt,et);

in = [pnWaves wavesfNames num2str(mon, '%02d') '.nc'];
outWaves = parseNCFile(in);
t = outWaves.var.time;
H = outWaves.var.waveHs;
Tp = outWaves.var.waveTp;
Hs = interp1(t,H,et);
fp = 1./interp1(t,Tp,et);
if isnan(Hs)
    Hs = 1;     % default values
end
if isnan(fp)
    fp = 1/8;
end

