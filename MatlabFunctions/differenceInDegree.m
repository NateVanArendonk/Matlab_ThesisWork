function degDiff = differenceInDegree(angle1,angle2)
% Calculates the difference between 2 normalized angles 
% https://stackoverflow.com/questions/32276369/calculating-absolute-differences-between-two-angles
% take from above link 

degDiff = -1*(normalizeDeg(angle1) - normalizeDeg(angle2));
end

function normDeg = normalizeDeg(angle)
% Normalize a angle to go from -180 to 180 
% https://stackoverflow.com/questions/32276369/calculating-absolute-differences-between-two-angles
% take from above link 
normDeg = @(x)(-mod(-angle+180,360)+180);
end

       