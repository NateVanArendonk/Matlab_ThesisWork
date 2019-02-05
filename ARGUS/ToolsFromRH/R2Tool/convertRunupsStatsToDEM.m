% read in all runups files from 18 Sept to end of Oct and find DEM values
% of extreme runup statistics.  These are the dates for which we have
% lidar data.

clear
pn = '/home/ruby/users/holman/research/RUNUP/discreteRunupTool/OUTPUTParkerFinal2/';
fns = dir([pn '*mat']);
pnOut = '/home/ruby/users/holman/research/RUNUP/discreteRunupTool/OUTPUTParkerFinal3/';

n = parseFilename(fns(1).name);                      % for later use
camInfo = findCamInfo(n.station, str2num(n.time));
xyzCam = camInfo(1).xyz;
yz0 = [850 0];
for i = 169: length(fns)
    i
    load([pn fns(i).name]);
    d = getDEMDay(runups.when);
    d.z = d.DEM;
    % find the plane slope around y = 850
    runups.RMax.xyzDEM = pixels2DEM(d.x, d.y, d.z, [runups.RMax.x runups.ym 0],xyzCam);
    runups.RMax.zWave = runups.RMax.xyzDEM(3) - runups.RMaxFromBrightest.env.zt;
    runups.R2.xyzDEM = pixels2DEM(d.x, d.y, d.z, [runups.R2.x runups.ym 0],xyzCam);
    runups.R2.zWave = runups.R2.xyzDEM(3) - runups.RMaxFromBrightest.env.zt;
    runups.R5.xyzDEM = pixels2DEM(d.x, d.y, d.z, [runups.R5.x runups.ym 0],xyzCam);
    runups.R5.zWave = runups.R5.xyzDEM(3) - runups.RMaxFromBrightest.env.zt;
    runups.R10.xyzDEM = pixels2DEM(d.x, d.y, d.z, [runups.R10.x runups.ym 0],xyzCam);
    runups.R10.zWave= runups.R10.xyzDEM(3) - runups.RMaxFromBrightest.env.zt;
    runups.R50.xyzDEM = pixels2DEM(d.x, d.y, d.z, [runups.R50.x runups.ym 0],xyzCam);
    runups.R50.zWave = runups.R50.xyzDEM(3) - runups.RMaxFromBrightest.env.zt;
    runups.RMaxFromBrightest.xyzDEM = pixels2DEM(d.x, d.y, d.z, [runups.RMaxFromBrightest.x runups.ym 0],xyzCam);
    runups.RMaxFromBrightest.zWave = runups.RMaxFromBrightest.xyzDEM(3) - runups.RMaxFromBrightest.env.zt;
    runups.RMax10Min.xyzDEM = pixels2DEM(d.x, d.y, d.z, [runups.RMax10Min.x runups.ym 0],xyzCam);
    runups.RMax10Min.zWave = runups.RMax10Min.xyzDEM(3) - runups.RMaxFromBrightest.env.zt;
    runups.RMaxFromBrightest.env.betaLidar = d.beta850;
    eval(['save ' pnOut fns(i).name ' runups'])
end
