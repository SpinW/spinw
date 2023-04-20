classdef mock_function < handle
    properties
       arguments = {};  % Arguments called with
       n_calls = 0;     % Number of times called
       func = '';       % Name of function
       filename = '';   % Name of function file
    end
    methods
        function mockobj = mock_function(function_name, return_value)
            if nargin < 2
                rv_str = '';
                return_value = '{[]}';
            else
                global mock_ret_val;
                if isempty(mock_ret_val)
                    mock_ret_val = struct();
                end
                if iscell(return_value)
                    mock_ret_val.(function_name) = return_value;
                else
                    mock_ret_val.(function_name) = {return_value};
                end
                rv_str = 'global mock_ret_val;';
                return_value = ['mock_ret_val.' function_name];
            end
            fnstr = [...
                     'function varargout = %s(varargin)\n' ...
                     '    persistent n_calls;\n' ...
                     '    persistent arguments;\n' ...
                     '    %s\n' ...
                     '    if nargin > 0 && ischar(varargin{1}) && strcmp(varargin{1}, ''check_calls'')\n' ...
                     '        varargout = {n_calls arguments};\n' ...
                     '        return;\n' ...
                     '    end\n' ...
                     '    if isempty(n_calls)\n' ...
                     '        n_calls = 1;\n' ...
                     '        arguments = {varargin};\n' ...
                     '    else\n' ...
                     '        n_calls = n_calls + 1;\n' ...
                     '        arguments = [arguments {varargin}];\n' ...
                     '    end\n' ...
                     '    if nargout > 0\n' ...
                     '        varargout = %s;\n' ...
                     '    end\n' ...
                     'end\n'];
            mockobj.func = function_name;
            mockobj.filename = sprintf('%s.m', function_name);
            fid = fopen(mockobj.filename, 'w');
            fprintf(fid, fnstr, function_name, rv_str, return_value);
            fclose(fid);
            whichfun = which(function_name);
            while ~strcmp(whichfun, fullfile(pwd, mockobj.filename))
                pause(0.1);
                whichfun = which(function_name);
            end
        end
        function delete(mockobj)
            delete(mockobj.filename);
            global mock_ret_val;
            if isfield(mock_ret_val, mockobj.func)
                mock_ret_val = rmfield(mock_ret_val, mockobj.func);
            end
        end
        function n_call = get.n_calls(mockobj)
            [n_call, ~] = feval(mockobj.func, 'check_calls');
            if isempty(n_call)
                n_call = 0;
            end
        end
        function arguments = get.arguments(mockobj)
            [~, arguments] = feval(mockobj.func, 'check_calls');
        end
    end
end
