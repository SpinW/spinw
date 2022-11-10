function has_toolbox = sw_hassymtoolbox()
% Checks if the running Matlab instance has the symbolic toolbox
%
% ### Syntax
%
% `has_toolbox = sw_hassymtoolbox()`
%
% ### Description
%
% This function checks if the symbolic toolbox is available
%

has_toolbox = license('test', 'symbolic_toolbox') == 1;
if has_toolbox
    try
        x = sym('x');
    catch ME
        if strcmp(ME.identifier, 'MATLAB:UndefinedFunction')
            has_toolbox = false
        else
            rethrow(ME)
        end
    end
end
