function T = XYZtoDist(x,y,z)

% Function that takes X,Y,Z points and calculates distance along a transect
% for them

%load file
swath.X = x;
swath.Y = y;
swath.Z = z;
%caluclate profile line
regswath = fitlm(swath.X, swath.Y);
mline = regswath.Coefficients.Estimate(2);
bline = regswath.Coefficients.Estimate(1);

%calculate point distance along profile
mperpline =  -1/mline;
Distance = zeros(length(swath.X), 1);
xprofpoint = zeros(length(swath.X), 1);
yprofpoint = zeros(length(swath.X), 1);
for i = 1:length(swath.X)
    
    bperpline = -(mperpline * swath.X(i))+ swath.Y(i);
    xprofpoint(i) = (bline - bperpline)/(mperpline - mline);
    yprofpoint(i) = (mperpline * xprofpoint(i)) + bperpline;
        
    pointmatr = [swath.X(1) swath.Y(1); xprofpoint(i) yprofpoint(i)];
    Distance(i) = pdist(pointmatr, 'euclidean');
end

%Create table of results
Elevation = swath.Z;
T = [Distance, Elevation];


