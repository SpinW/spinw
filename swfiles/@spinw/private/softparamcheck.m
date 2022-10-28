function softparamcheck(params_to_check, func_name, param, varargin)
% Checks if any parameters have been provided in varargin, but set to
% empty, and raises an error if so. This function is needed because by
% default 'soft' params will silently be set to empty if their shape is
% incorrect, causing unexpected behaviour for users.
%
% Input:
%
% params_to_check    A list of strings of the parameters to check e.g.
%                    ["S", "k"]
% func_name          The name of the function this is being called from
%                    to be used in the error identifier text
% param              The output of sw_readparam. If an argument exists in
%                    varargin, but has been set to empty in param, we know
%                    it has been silently ignored so raise an error
% varargin           Varargin that was used as input to sw_readparam to
%                    create param, this should be name-value pairs or a
%                    struct
%
if isempty(varargin)
    return
end
names = vararginnames(varargin{:});
err_str = [];
for i = 1:length(params_to_check)
    if any(strcmpi(names, params_to_check(i))) ...
       && isempty(param.(params_to_check(i)))
        err_str = [err_str params_to_check(i)];
    end
end
if ~isempty(err_str) > 0
    error(['spinw:' char(func_name) ':WrongInput'], ...
          'Incorrect input size for ' + join(err_str, ', '));
end