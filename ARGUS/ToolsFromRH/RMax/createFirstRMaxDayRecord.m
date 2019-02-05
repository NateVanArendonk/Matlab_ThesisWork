function RMaxDay = createFirstRMaxDayRecord(fn, camInfo)
%   RMaxDay = createFirstRMaxDayRecord(fn, camInfo)
%
% Given a filename of a single brightest image, create a mostly empty
% RMaxDay structure and manually create a seed shoreline for subsequent
% analyses.  This routine would be called to initiate a sequence of
% analyses.  Thereafter, each new RMaxDay is built on the previous day so
% this routine is no longer needed.
% fn can be any representative brightest image (including the one
% that will be used for the first analysis).  It should be appropriate for
% digitizing a seed shoreline, i.e. the runup edge in the brightest images
% should be easy to find.  The shoreline seed is created manually, by
% clicking approximate shoreline locations in each of a set of brightest
% images, one for each camera.  The seed is also saved locally as
% seedShoreline.mat simply to avoid having to re-digitize during
% development.  
% camInfo contains the list of cam numbers, their xyz locations, their
% geometries as m-vectors and their cam structures

% get camera number list and metadata
n = parseFilename(fn);
eval([n.station 'BrightestParams'])     % may need to be better implemented
et = str2num(n.time);
camNums = [camInfo.num];
TZOffset = DBGetStationTZOffset(n.station);  % in minutes to standard time

% create empty struct
RMaxDay = createEmptyRMaxDayStruct(length(camNums), params);
RMaxDay.camNums = camNums;
foo = [camInfo.xyz];
RMaxDay.xyzCams = reshape(foo',3,[])';

% get environmental info for temporary use to find x0 shore
eval(['[zt, Hs, fp] = ' params.envInputProg '(et);']);
eval(['R2 = ' params.findZEstProg '(Hs,fp,params.beta);'])
delx = (zt+R2)/params.beta;

% get list of files that are synchronous with fn
for i = 1: length(camNums)
    n.camera = camNums(i);
    fnList(i,:) = findArgusImages(n);
end

% Look for seedShoreline that was previously digitized.  If none or if too
% old, digitize again.
maxGapInDays = 14;      % digitize new seed shoreline if larger time gap 
if exist('seedShoreline.mat')
    load('seedShoreline')
    dt = seed.when - str2num(fn(1:10));
    if abs(dt) > maxGapInDays*24*3600;
        clear seed
    end
end
if ~exist('seed')
    xyz0 = digitizeBrightestShoreline(fnList, delx, params, camInfo);
    seed.xyzSeed = xyz0;
    seed.when = str2num(fnList(1,1:10));
    save seedShoreline seed
else
    load('seedShoreline')
end

% now add the new elements to the structure RMaxDay
Ny = length(params.ys);
RMaxDay.whenPrevious = et - 24*3600;  % fake previous day
RMaxDay.seed.Kalman.xsKal = seed.xyzSeed(:,1);
RMaxDay.seed.Kalman.PNew = repmat(params.xsSig.^2, Ny, 1);
RMaxDay.seed.Kalman.when = et - (24*3600);  % back date one day
RMaxDay.beta.Kalman.betaKal = repmat(params.beta,Ny,1);
RMaxDay.beta.Kalman.PNew = repmat(params.betaSig^2, Ny, 1);
RMaxDay.mapped.zEst.Kalman.when = et - 24*3600;    % fake time
