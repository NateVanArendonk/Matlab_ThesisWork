function camInfo = findCamInfo(stationNameStr, epoch)
%   camInfo = findCamInfo(stationNameStr, epoch)
%
% Find all of the important camera info for a site at an epoch.  For the
% CIL this is done by database calls, others will have to provide their own
% info.  Required info is returned in the structure camInfo containing
% fields
%   camNums         - NCam by 1 sorted list of camera numbers
%   xyzCams         - NCam by 3 table of xyz locations
%   m               - NCam by 11 table of geometries
%   camStruct       - database return of camera struct with distortion info

a = DBGetCamerasByStation(stationNameStr, epoch);
camNumber = [a.cameraNumber];
xyzCam = [[a.x]' [a.y]' [a.z]'];
[camNumber, ind] = sort(camNumber);
xyzCam = xyzCam(ind,:);
for i = 1: length(camNumber)
    camInfo(i).num = camNumber(i);
    camInfo(i).xyz = xyzCam(i,:);
    pick = find([a.cameraNumber]==camNumber(i));
    foo = DBGetCurrentGeom(a(pick).id, epoch);
    camInfo(i).m = foo.m;
    camInfo(i).cam = a(pick);
    camInfo(i).ip = DBLegacyIP(a(pick));
end