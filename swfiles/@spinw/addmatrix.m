function addmatrix(obj, varargin)
% adds new 3x3 matrix
% 
% ### Syntax
% 
% `addmatrix(obj,Name,Value)`
% 
% ### Description
% 
% `addmatrix(obj,Name,Value)` adds a new $[3\times 3]$ matrix to the
% [spinw.matrix] field of `obj`. The added matrices can be later assigned
% to bonds, single ion anisotropy terms or g-tensors of magnetic atoms. If
% the given matrix label already exists in `obj`, instead of adding new
% matrix the existing one will be overwritten.
% 
% ### Examples
% 
% The first example adds a diagonal matrix `eye(3)`, that can describe
% Heisenberg interaction if assigned to a bond. The second example adds an
% ansisymmetric matrix that can decribe Dzyaloshinskii-Moriya (DM)
% interaction if assigned to a bond.
%
% ```
% >>crystal = spinw
% >>crystal.addmatrix('value',1,'label','J_1')
% >>crystal.matrix.mat>>
% >>crystal.addmatrix('value',[1 0 0],'label','J_1')
% >>crystal.matrix.mat>>
% ```
% 
% ### Input Arguments
% 
% `obj`
% : [spinw] object.
% 
% ### Name-Value Pair Arguments
% 
% `'value'`
% : The actual numerical values to be added as a matrix. It can have the
%   following shapes:
%   * $[3\times 3]$ the given values will be stored in [spinw.matrix] as
%     they are given.
%   * $[1\times 1]$ the given value will be multiplied with `eye(3)`.
%   * `[Mx My Mz]` the given triplet will be used to define an
%     antisymmetric matrix `M = [0 M3 -M2;-M3 0 M1;M2 -M1 0]`. 
% 
% `'label'`
% : Label string for plotting default value is `'matI'`, where $I$ is the index
%   of the matrix. Add '-' to the end of the label to plot bond as dashed
%   line/cylinder.
% 
% `'color'`
% : Color for plotting, either row vector
%   that contains color RGB codes (values of 0-255), or a string with the
%   name of the color, for possible colors names [swplot.color]. Default
%   color is a random color.
% 
% ### Output Arguments
% 
% The `obj` output will contain the additional matrix in the [spinw.matrix]
% field.
% 
% ### See Also
% 
% [spinw] \| [swplot.color]
%
% *[DM]: Dzyaloshinski-Moriya
% *[RGB]: Red-Green-Blue
%

if nargin < 2
    swhelp spinw.addmatrix
    return
end

inpForm.fname  = {'value'    'mat'      'label' 'color' };
inpForm.defval = {[]         []         {}      []      };
inpForm.size   = {[-1 -2 -3] [-4 -5 -6] [-7 -8] [-9 -10] };
inpForm.soft   = {true       true       true    true    };

newMat = sw_readparam(inpForm, varargin{:});

if ~isnumeric(newMat.value) && ~isa(newMat.value,'sym')
    error('spinw:addmatrix:WrongInput','Matrix value has to be numeric or symbolic variable!')
end

if ~isempty(newMat.value)
    newMat.mat = newMat.value;
end

if numel(newMat.mat) == 1
    newMat.mat = newMat.mat*eye(3);
end

if numel(newMat.mat) == 3
    M = newMat.mat;
    newMat.mat = [0 M(3) -M(2);-M(3) 0 M(1);M(2) -M(1) 0];
end


% Defult Heisenberg matrix.
if isempty(newMat(1).mat)
    if ~iscell(newMat.label)
        nLabel = 1;
    else
        nLabel = size(newMat.label,2);
    end
    nJ = nLabel;
    if nJ == 0
        error('spinw:addmatrix:WrongInput','Define some matrix property!');
    else
        newMat.mat = repmat(eye(3),[1 1 nJ]);
        if ~obj.symbolic
            warning('spinw:addmatrix:NoValue','No valid value was given for the new matrix, default value used!');
        end
    end
end

if ~isfield(newMat,'label') || isempty([newMat.label])
    idx = size(obj.matrix.mat,3)+1;
    for ii = 1:numel(newMat)
        newMat(ii).label = {};
        for jj = 1:size(newMat(ii).mat,3)
            newMat(ii).label = [newMat(ii).label {['mat' num2str(idx)]}];
            idx = idx + 1;
        end
    end
end

if ~isfield(newMat,'color') || isempty([newMat.color])
    for ii = 1:numel(newMat)
        %newMat(ii).color = repmat([255;0;0],1,size(newMat(ii).mat,3));
        % generate random color
        newMat(ii).color = repmat(swplot.color(randi(141),1),1,size(newMat(ii).mat,3));
    end
else
    for ii = 1:numel(newMat)
        newMat(ii).color = swplot.color(newMat(ii).color);
    end
    
end

% If the structure is a vector concatenates it.
newJVect     = newMat;
newMat       = struct;
newMat.mat   = [];
newMat.label = {};
newMat.color = [];

for ii = 1:numel(newJVect)
    %if ~any(size(newJVect(ii).color)-[1 3])
    %    newJVect(ii).color = newJVect(ii).color';
    %end
    if ~iscell(newJVect(ii).label)
        newJVect(ii).label = {newJVect(ii).label};
    end
    
    newMat.mat   = cat(3,newMat.mat,newJVect(ii).mat);
    newMat.label = [newMat.label newJVect(ii).label];
    newMat.color = [newMat.color newJVect(ii).color];
end

% Goes through the elements one-by-one.
for ii = 1:size(newMat.mat,3)
    newJItem.mat   = newMat.mat(:,:,ii);
    newJItem.label = newMat.label(ii);
    newJItem.color = int32(newMat.color(:,ii));
    newObj.matrix  = newJItem;
    
    %spinw.validate(newObj,'matrix');
    
    cIdx = find(strcmp(obj.matrix.label,newJItem.label));
    
    if numel(cIdx)>1
        error('spinw:addmatrix:LabelError','There are multiple matrices with the same label!');
        
    elseif numel(cIdx) == 1
        % Replace the coupling type with the new one.
        if obj.symb && ~isa(newJItem.mat,'sym')
            obj.matrix.mat(:,:,cIdx) = newJItem.mat*sym(newJItem.label,'real');
        else
            obj.matrix.mat(:,:,cIdx) = newJItem.mat;
        end
        obj.matrix.label(cIdx)   = newJItem.label;
        obj.matrix.color(:,cIdx) = newJItem.color;
        
    else
        % Add the coupling type to the list.
        if (obj.symb) && ~isa(newObj.matrix.mat,'sym')
            for jj = 1:size(newObj.matrix.mat,3)
                if newObj.matrix.label{jj}(end) == '-'
                    symVar = sym(newObj.matrix.label{jj}(1:(end-1)),'real');
                else
                    symVar = sym(newObj.matrix.label{jj},'real');
                end
                obj.matrix.mat = cat(3,obj.matrix.mat,newObj.matrix.mat(:,:,jj)*symVar);
            end
        else
            obj.matrix.mat = cat(3,obj.matrix.mat,newObj.matrix.mat);
        end
        obj.matrix.label = [obj.matrix.label newObj.matrix.label];
        obj.matrix.color = [obj.matrix.color newObj.matrix.color];
        
    end
end

%spinw.validate(obj);

end
