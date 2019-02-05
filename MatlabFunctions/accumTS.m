%% accumTS- aggregate data to users specified columns of data
%
% Function uses matlab function ACCUMARRAY to accumulate timeseries data to
% user specified input data of uniqueness.
%
% Expected inputs is a timeseries of data with all but the last column used
% for generate unique cominations, and the last column is where the
% summariziation is done.
%
% datain = r x c
% fn = supplied string specifying a valid summary statistic function.
% Example: @sum, @max, @min, etc.
% 
% dataout = r x c
% 
% Written by:
%   Jeff Burkey
%   King County- DNRP
%   email: jeff.burkey@kingcounty.gov
%   June 9, 2009
% 
% Syntax:
%   [dataout] = accumTS(datain,fn)
%     d = accumTS(din,@max);
function [dataout] = accumTS(datain,fn)
    % Aggregate values using the user specified function.
    [b1, m1, n1] = unique(datain(:,1:end-1),'rows');
    dataout = [b1 accumarray(n1,datain(:,end),[],fn)];
end