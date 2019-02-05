% RMax Toolbox for detection of maximum runup from brightest images
%
% User Routines
%   testBrightestExtract    - Master routine
%   argus02bBrightestParams - analysis parameters for argus02b (example)
%   argus02bQCLimits        - quality control ranges for argus02b (example)
%
% Major Functions
%   createFirstRMaxRecord   - initialize structure and find seed shoreline
%   digitizeBrightestShoreline - manually digitize seed frame
%   findRMax                - analyze set of images for one instant in time
%   makeMapsRMaxDay         - make foreshore map and new seed for one day
%
% Debug/Visualization 
%   examineRMaxDay          - plots all points, interp update and Kalman
%   showRMaxUpDateResults   - scatter plot of days's update results
%   showRMaxKalmanResults   - maps of update and Kalman results
%   showRMaxSeedAndBetaResults - shows shoreline and beta results for a day
%   
% Support Functions
%   createEmptyRMaxDayStruct - create an empty RMaxDay structure
%   getDuckEnvirInfo        - find tide, Hs, fp for collection (example)
%   findCamInfo             - return camInfo structure
%   StockdonEquation        - estimate R2 from Stockdon
%   findBrightEdge          - RMax algorithm engine
%   brightestQC             - quality control on estimates
%   findxyzDEM              - map RMax locations onto a DEM
%   findxyzEst              - map RMax locations onto an estimated z-level
%   findClosestDEM          - retrieve a DEM pathname for any time
%   getLowTideDEM           - get DEM at nearest low tide
%   pixels2DEM              - map pixels locations onto a DEM
%   mergeResults            - combine estimates from different cameras
%   sigmoid1                - 1D sigmoid function for fitting
%   findQ                   - find process error, for any wave height
%   KalmanFilter            - do Kalman filtering
%   getRMaxVersion          - version number source

