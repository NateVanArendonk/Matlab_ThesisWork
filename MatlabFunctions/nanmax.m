function [varargout] = nanmax(varargin)
% Minimum value, ignoring NaNs


[varargout{1:nargout}] = max(varargin{:});
end