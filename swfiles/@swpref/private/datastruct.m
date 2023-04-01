function d = datastruct()
% Function called to create a default preference object.
%
%  {{warning Internal function for the Spin preferences.}}
%
% ### Syntax
%
% 'prefs = datastruct()'
%
% ### Description
%
% 'prefs = datastruct()' creates a structure with the following fields:
%
% 'Name' a cell array giving the name of each dynamic property.
%
% 'Validation' a cell array with functions to evaluate when a
% property is set.
%
% 'Value' a cell array giving the default value of a dynamic property
%
% 'Label' a cell array giving a short description on the dynamic
% property.
%

d.Name = {
    'fid',...
    'expert',...
    'tag',...
    'nmesh',...
    'maxmesh',...
    'npatch',...
    'fontsize',...
    'tid',...
    'colormap',...
    'usemex',...
    'docurl'
    };

Size = {
    [1, 1],...
    [1, 1],...
    [    ],...
    [1, 1],...
    [1, 1],...
    [1, 1],...
    [1, 1],...
    [1, 1],...
    [1, 1],...
    [1, 1],...
    [    ]
    };

d.Validation = {
    {@isnumeric, @mustBeInteger, @mustBeNonnegative, @(x) check_size(x,Size{1}), @(x) mustBeLessThan(x,256)},...
    {@islogical, @(x) check_size(x,Size{2})},...
    {@ischar},...
    {@isnumeric, @mustBeInteger, @mustBeNonnegative, @(x) check_size(x,Size{4}), @(x) mustBeGreaterThanOrEqual(x,1), @(x) mustBeLessThan(x,256)},...
    {@isnumeric, @mustBeInteger, @mustBeNonnegative, @(x) check_size(x,Size{5}), @(x) mustBeGreaterThan(x,1), @(x) mustBeLessThan(x,256)},...
    {@isnumeric, @mustBeInteger, @mustBeNonnegative, @(x) check_size(x,Size{6}), @(x) mustBeGreaterThan(x,1), @(x) mustBeLessThan(x,256)},...
    {@isnumeric, @mustBeInteger, @mustBeNonnegative, @(x) check_size(x,Size{7}), @(x) mustBeGreaterThan(x,4), @(x) mustBeLessThan(x,256)},...
    {@isnumeric, @mustBeInteger, @mustBeNonnegative, @(x) check_size(x,Size{8}), @(x) mustBeLessThan(x,256)},...
    {@(x) check_size(x,Size{9}), @(x) isa(x,'function_handle')},...
    {@(x) check_size(x,Size{10}), @(x) check_mex(x), @islogical},...
    {@ischar, @(x) strfind(x,'http')}
    };

d.Value = {
    1,...
    false,...
    'swplot',...
    1,...
    6,...
    20,...
    12,...
    1,...
    @cm_inferno,...
    false,...
    'https://spinw.github.io/spinwdoc'
    };

d.Label =  {
    'file identifier for text output, default value is 1 (Command Window)'...
    'expert mode (1) gives less warnings (not recommended), default value is 0'...
    'defines the tag property of the crystal structure plot figures'...
    'default number of subdivision of the icosahedron for plotting'...
    'maximum number of subdivision of the icosahedron for plotting'...
    'number of edges for patch object'...
    'fontsize for plotting'...
    'identifier how the timer is shown, default value is 1 (Command Window), value 2 gives graphical output'...
    'default colormap'...
    'if true, mex files are used in the spin wave calculation'...
    'url to the documentation server'...
    };

    function out = check_size(obj,S)
        % checks to see if an object is the wrong size.
        %
        %  {{warning Internal function for the Spin preferences.}}
        %
        % ### Syntax
        %
        % 'logical = check_size(toBeChecked,size)'
        %
        % ### Description
        %
        % 'logical = check_size(toBeChecked,size)' checks to see if an
        % object 'obj 'is the expected size given by 'size'. An error is
        % thrown if there is a difference.
        %
        
        sz = size(obj);
        if ~all(sz == S)
            error('spref:WrongSize','Value to be asigned is the wrong size [%i, %i] not [%i, %i]',sz(1), sz(2), S(1), S(2))
        else
            out = 1;
        end
    end

    function out = check_mex(val)
        % checks to see if mex files are available.
        %
        %  {{warning Internal function for the Spin preferences.}}
        %
        % ### Syntax
        %
        % 'logical = check_mex(obj)'
        %
        % ### Description
        %
        % 'logical = check_mex(obj)' checks to see if files 'chol_omp' and
        % 'eig_omp' are present in the MATLAB path.An error is thrown if
        % they do not exist.
        %
        
        % Do the checks only if we are trying to set usemex = true
        if val == 0
            out = 0;
            return
        end
        
        if ~(exist('chol_omp','file')==3 && exist('eig_omp','file')==3)
            % There is a path error for < R2017a
            if (exist('chol_omp','file')==2 && exist('eig_omp','file')==2)
                p = which('chol_omp');
                if ~(exist('chol_omp','file')==3 && exist('eig_omp','file')==3)
                    error('spref:MissingMex','Necessary mex files are missing, compile them!')
                end
            else
                error('spref:MissingMex','Necessary mex files are missing, compile them!')
            end
        else
            out = 1;
        end
    end

    function mustBeInteger(A)
        % Validate that value is integer or issue error
        %
        %  {{warning Internal function for the Spin preferences.}}
        %
        %   MUSTBEINTEGER(A) issues an error if A contains non integer values.
        %   A value is integer if it is real, finite, and equal to the result
        %   of taking the floor of the value.
        %
        %   Modified from mathworks verion for < R2017a
        
        ME = MException('MATLAB:validators:mustBeNumericOrLogical','Value must be integer.');
        if ~all(isnumeric(A) || islogical(A))
            throwAsCaller(ME)
        end
        
        if ~isreal(A)
            throwAsCaller(ME)
        end
        
        if ~all(isfinite(A(:))) || ~all(A(:) == floor(A(:)))
            throwAsCaller(ME)
        end
    end

    function mustBeNonnegative(A)
        % Validate that value is nonnegative or issue error
        %
        %  {{warning Internal function for the Spin preferences.}}
        %
        %   MUSTBENONNEGATIVE(A) issues an error if A contains negaitive values.
        %   A value is nonnegative if it is greater than or equal to zero.
        %
        %   Modified from mathworks verion for < R2017a
        
        if ~all(A(:) >= 0)
            ME = MException('MATLAB:validators:mustBeNonnegative','Value/s must be positive.');
            throwAsCaller(ME)
        end
    end
    function mustBeGreaterThan(A, B)
        % Validate that value is greater than a specified value or issue error
        %
        %  {{warning Internal function for the Spin preferences.}}
        %
        %   MUSTBEGREATERTHAN(A,B) issues an error if A is not greater than B.
        %   MATLAB calls gt to determine if A is greater than B.
        %
        %   Modified from mathworks verion for < R2017a
        
        if ~all(A(:) > B)
            ME = MException('MATLAB:validators:mustBeGreaterThan','All values must be greater than %f.',B);
            throwAsCaller(ME)
        end
    end
    function mustBeLessThan(A, B)
        % Validate that value is less than a specified value or issue error
        %
        %  {{warning Internal function for the Spin preferences.}}
        %
        %   MUSTBELESSTHAN(A,B) issues an error if A is not less than B.
        %   MATLAB calls gt to determine if A is less than B.
        %
        %   Modified from mathworks verion for < R2017a
        
        if ~all(A(:) < B)
            ME = MException('MATLAB:validators:mustBeLessThan','All values must be less than %f.',B);
            throwAsCaller(ME)
        end
    end

    function mustBeGreaterThanOrEqual(A, B)
        % Validate that value is greater than or equal to a specified value or issue error
        %
        %  {{warning Internal function for the Spin preferences.}}
        %
        %   MUSTBEGREATERTHANOREQUAL(A,B) issues an error if A is not greater than or equal to B.
        %   MATLAB calls gt to determine if A is greater than or equal to B.
        %
        %   Modified from mathworks verion for < R2017a
        
        if ~all(A(:) >= B)
            ME = MException('MATLAB:validators:mustBeGreaterThanOrEqual','All values must be greater than or equal to %f.',B);
            throwAsCaller(ME)
        end
    end

    function mustBeLessThanOrEqual(A, B)
        % Validate that value is less than or equal to a specified value or issue error
        %
        %  {{warning Internal function for the Spin preferences.}}
        %
        %   MUSTBELESSTHANOREQUAL(A,B) issues an error if A is not less or equal than B.
        %   MATLAB calls gt to determine if A is less than oir equal B.
        %
        %   Modified from mathworks verion for < R2017a
        
        if ~all(A(:) <= B)
            ME = MException('MATLAB:validators:mustBeLessThanOrEqual','All values must be less than or equal to %f.',B);
            throwAsCaller(ME)
        end
    end
end
