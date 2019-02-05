function [varargout] = nanmin(varargin)
% Minimum value, ignoring NaNs


[varargout{1:nargout}] = min(varargin{:});
end