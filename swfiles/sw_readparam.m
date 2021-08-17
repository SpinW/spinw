function [input, remainder] = sw_readparam(format, varargin)
% parse input arguments
%
% ### Syntax
%
% `parsed = sw_readparam(format,Name,Value)`
%
% ### Description
%
% `parsed = sw_readparam(format,Name,Vale)` parses name-value pair
% arguments. The parsing is controlled by the given `format` input. The
% name-value pairs are converted into the parsed struct which has field
% names identical to the given parameter names and corresponding values
% taken from the input. `format` can also define required dimensionality of
% a given value and default values for select parameters.
%
% `sw_readparam` is used in most of the method functions of [spinw].
%
% ### Input Arguments
%
% `format`
% : A struct with the following fields:
%   * `fname` Field names, $n_{param}$ strings in cell.
%   * `size` Required dimensions of the corresponding value in a cell of
%     $n_{param}$ vectors. Negative integer means dimension has to match with
%     any other dimension which has the identical negative integer.
%   * `defval` Cell of $n_{param}$ values, provides default values for
%     missing parameters.
%   * `soft` Cell of $n_{param}$ logical values, optional. If `soft(i)` is
%     true, in case of missing parameter value $i$, no warning will be
%     given.
%

if nargin == 0
    swhelp sw_readparam
    return
end

if nargin == 2 && ischar(varargin{1}) && varargin{1} == '-'
    % return all variable strings onto the Command Window
    varStr = sprintf('%s ',format.fname{:});
    clipboard('copy',varStr(1:end-1));
    warning('sw_readparam:ClipboardCopy','Input parameter names are copied to the clipboard!');
end

% create a showWarn field to check whether to show warnings (default true)
%format.fname  = [format.fname  {'showWarn'}];
%format.defval = [format.defval {true      }];
%format.size   = [format.size   {[1 1]     }];
%if isfield(format,'soft')
%    format.soft   = [format.soft   {true      }];
%end

check = true;
if nargin==3 && isstruct(varargin{1}) && strcmp(varargin{2},'noCheck')
    % just return the structure without checking for speedup, use it only
    % for internal calls
    raw   = varargin{1};
    check = false;
elseif (nargin>2) && (mod(nargin,2) == 1)
    nPar = nargin-1;
    raw = struct;
    for ii = 1:2:nPar
        raw.(varargin{ii}) = varargin{ii+1};
    end
elseif nargin == 2
    raw = varargin{1};
elseif nargin == 1
    raw = struct;
else
    MException('sw_readparam:WrongParameter','Parameter name-value pairs are expected!').throwAsCaller;
end

if ~isstruct(raw)
    MException('sw_readparam:WrongParameter','Parameter name-value pairs are expected!').throwAsCaller;
end

fName     = format.fname;
rName     = fieldnames(raw);
storeSize = zeros(20,1);
input     = struct;

usedField = false(1,numel(rName));

% Go through all fields.
for ii = 1:length(fName)
    
    rawIdx = find(strcmpi(rName,fName{ii}));
    
    if any(rawIdx)
        rawIdx = rawIdx(1);
        usedField(rawIdx) = true;
        
        inputValid = true;
        
        % Go through all dimension of the selected field to check size.
        if check
            for jj = 1:length(format.size{ii})
                if format.size{ii}(jj)>0
                    if format.size{ii}(jj) ~= size(raw.(rName{rawIdx}),jj)
                        inputValid = false;
                    end
                else
                    if storeSize(-format.size{ii}(jj)) == 0
                        storeSize(-format.size{ii}(jj)) = size(raw.(rName{rawIdx}),jj);
                    else
                        if storeSize(-format.size{ii}(jj)) ~= size(raw.(rName{rawIdx}),jj)
                            inputValid = false;
                        end
                        
                    end
                end
            end
        end
        if inputValid
            input.(fName{ii}) = raw.(rName{rawIdx});
        else
            if isfield(format,'soft') && format.soft{ii}
                input.(fName{ii}) = format.defval{ii};
            else
                MException('sw_readparam:ParameterSizeMismatch',['Input parameter size mismatch in parameter ''' fName{ii} '''!']).throwAsCaller;
            end
        end
    else
        if isfield(format,'defval') && (any(size(format.defval{ii})) || (isfield(format,'soft') && format.soft{ii}))
            input.(fName{ii}) = format.defval{ii};
        else
            MException('sw_readparam:MissingParameter',['Necessary input parameter ''' fName{ii} ''' is missing!']).throwAsCaller;
        end
    end
end

if ~all(usedField)
    mName = rName(~usedField);
    if nargout > 1
        remainder = {};
        for ii = 1:numel(mName)
            remainder =  [remainder {mName{ii} raw.(mName{ii})}];
        end
    else
        if numel(mName) == 1
            wName = [' "' mName{1} '"'];
        elseif numel(mName) == 2
            wName = ['s "' mName{1} '" and "' mName{2} '"'];
        else
            wName = ['s ' sprintf('"%s", "',mName{1:(end-2)}) mName{end-1} '" and "' mName{end} '"'];
        end
        warning('sw_readparam:UnreadInput','Unregognised input parameter%s!',wName);
    end
elseif nargout > 1
    remainder = {};
end

end
