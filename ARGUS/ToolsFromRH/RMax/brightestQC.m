function QC = brightestQC(results, QCLims)
%   QC = brightestQC(results, QCLims);
%
% returns a vector of flags for acceptable brightest estimates.
% Criteria are contained in the matrix QCLims which describes the lower and
% upper limits for each of the eight variables in results

for i = [1:size(results,2)-1]
    good(:, i) = ((results(:,i) > QCLims(i,1)) & (results(:,i) < QCLims(i,2)));
end
good(:,4) = 1;      % force camera number to not be limiting yet
QC = all(good')';
QC(isnan(QC)) = 0;

