function y = nansum(varargin)
% sum of data, ignorning nan's

narginchk(1,2);
y = sum(varargin{:},'omitnan');
end