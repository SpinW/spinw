function names = vararginnames(varargin)
% Given varargin, returns a list of the given
% argument names (so we know which arguments a user has passed to a function).
% Note that inputs to SpinW functions can either be name-value pairs or a
% struct
%
% Input:
%
% varargin    Variable-length argument list (name value pairs) or struct
if isstruct(varargin{1})
    varargin_struct = varargin{1};
else
    varargin_struct = struct(varargin{:});
end
names = fields(varargin_struct);
end