function addcoupling(obj, varargin)
% assigns an exchange matrix to a bond
% 
% ### Syntax
% 
% `addcoupling(obj,Name,Value)`
% 
% ### Description
% 
% `addcoupling(obj,Name,Value)` assigns a matrix (will be used as exchange
% matrix) to a given bond after bonds are generated using
% [spinw.gencoupling].
% 
% ### Examples
% 
% To add the $J_1$ diagonal matrix to all second neighbor bonds
% between magnetic atoms use the following:
%
% ```
% >>cryst = sw_model('squareAF',1)
% >>cryst.addmatrix('label','J1','value',diag([1 0.1 0.1]))
% >>cryst.gencoupling
% >>cryst.addcoupling('mat','J1','bond',2)
% >>plot(cryst,'range',[2 2 1])
% >>snapnow
% ```
% 
% ### Input Arguments
% 
% `obj`
% : [spinw] object.
% 
% ### Name-Value Pair Arguments
% 
% `'mat'`
% : Label (string) or index (integer) of the matrix that will be assigned to
%   selected bonds, e.g. `'J1'`.
% 
% `'bond'`
% : Integer that selects bonds, e.g. 1 for first neighbor, 2 for second
%   neighbor, etc. The given value is compared to the `obj.coupling.idx`
%   vector and the exchange matrix will be assigned to matching bonds.
%   `'bond'` can be also a row vector to assign matrices to multiple bonds.
% 
% `'atom'`
% : Contains labels of atoms (string) or index of atoms (integer) that is
%   compared to [spinw.unit_cell] where all symmetry inequivalent atoms are
%   stored. If a single string label or number is given, e.g. `'Cr1'` only
%   Cr1-Cr1 bonds will be assigned. If a cell with 2 strings, e.g. `{'Cr1'
%   'Cr2'}` only Cr1-Cr2 bonds will be assigned. Default value is `[]`.
% 
% `'subIdx'`
% : If the above options are not enough to select the desired
%   bonds, using `subIdx` bonds can be selected one-by-one from
%   the list of bonds that fulfill the constraint of `atom` and `bond`.
% 
% `'type'`
% : Type of the coupling with possible values of:
%   * `'quadratic'`     Quadratic exchange, default.
%   * `'biquadratic'`   Biquadratic exchange.
% 
% `'sym'`
% : If `true`, symmetry operators will be applied on the exchange
%   matrices to generate the coupling on symmetry equivalent
%   bonds, if `false` all symmetry equivalent bonds will have the
%   same exhcange matrix.
% 
% {{warning Setting `atom` or `subIdx` parameters will remove the symmetry
% operations on the selected bonds. This means that assigning any
% non-Heisenberg exchange matrix will break the space group defined in
% `obj.lattice.sym`. Effectively reducing the symmetry of the given bond to
% `P0`}}
%
% ### Output Arguments
% 
% The function adds extra entries to the [spinw.coupling] property of
% `obj`. Specifically it will modify `obj.coupling.mat_idx`,
% `obj.coupling.type` and `obj.coupling.sym` matrices.
% 
% ### See Also
% 
% [spinw] \| [spinw.gencoupling] \| [spinw.addmatrix]


if isempty(obj.matom.r)
    error('spinw:addcoupling:NoMagAtom','There is no magnetic atom in the unit cell with S>0!');
end

inpForm.fname  = {'mat'  'bond' 'atom' 'subIdx' 'type' 'sym' };
inpForm.defval = {[]     []     []      []      []     []    };
inpForm.size   = {[1 -1] [1 -2] [-3 -4] [1 -5]  [1 -6] [1 -7]};
inpForm.soft   = {false  false  true    true    true   true  };

param = sw_readparam(inpForm, varargin{:});

if ~isnumeric(param.mat)
    param.mat = find(ismember(obj.matrix.label,param.mat));
end

if isempty(param.type)
    param.type = 0*param.mat;
end

if ischar(param.type)
    param.type = {param.type};
end

if iscell(param.type)
    % check string input for the coupling type
    qSel  = strcmp(param.type,'quadratic');
    bqSel = strcmp(param.type,'biquadratic');
    param.type = nan(size(param.type));
    param.type(qSel)  = 0;
    param.type(bqSel) = 1;
    if any(isnan(param.type))
        error('spinw:addcoupling:WrongInput',['Wrong coupling type, '...
            'currently only quadratic and biquadratic exchanges '...
            'are supported!']);
    end
    
end

if isempty(param.sym)
    nosympar = true;
    if isempty(param.subIdx)
        param.sym = true;
    else
        param.sym = false;
    end
else
    nosympar = false;
end

if numel(param.sym) == 1
    param.sym = param.mat*0+param.sym;
end

if any(size(param.mat)~=size(param.type))
    error('spinw:addcoupling:WrongInput',['A coupling type has to be '...
        'provided for each input matrix!'])
end

if ~isempty(param.atom)
    % select atoms based on label or index
    if ~isnumeric(param.atom)
        % atom labels provided - convert to index
        if ~iscell(param.atom)
            % single label provided as string - convert to cell for consistency
            param.atom = {param.atom};
        end
        % replace label with index (note can have more than one atom with
        % the same label)
        param.atom = cellfun(@(label) find(strcmp(obj.unit_cell.label, label)), ...
            param.atom, 'UniformOutput', false); % outputs cell
        if any(cellfun(@(idx) isempty(idx), param.atom))
            error('spinw:addcoupling:WrongInput', 'Atom label does not exist in unit cell.');
        end
    else
        % indices provided - convert to a cell if not already
        if ~iscell(param.atom)
            param.atom = num2cell(param.atom);
        end
        % check that all are present in matom
        if ~all(cellfun(@(idx) any(obj.matom.idx==idx), param.atom))
            error('spinw:addcoupling:WrongInput', 'Atom index does not correspond to a valid magentic atom.');
        end
    end
    if numel(param.atom)>2
        error('spinw:addcoupling:WrongInput','A maximum of 2 atom labels can be provided.');
    else
        aIdx1 = param.atom{1}; % atom 1 index in matom.idx
        if numel(param.atom)>1
           aIdx2 = param.atom{2}; % atom 2 index in matom.idx
        else
           aIdx2 = aIdx1;
        end
    end
end

if ~isempty(param.atom) || ~isempty(param.subIdx)
    if numel(param.bond) > 1
        warning('spinw:addcoupling:CouplingSize',['bond parameter has to be '...
            'scalar, only the first bond is selected!']);
        param.bond = param.bond(1);
    end
    %     if obj.sym
    %         error('spinw:addcoupling:SymmetryProblem',['atom and subIdx options are not allowed '...
    %             'when crystal symmetry is not P0!']);
    %     end
end

idx = ismember(obj.coupling.idx,param.bond);
if ~any(idx)
    error('spinw:addcoupling:CouplingError',['Coupling with idx=%d does '...
        'not exist, use gencoupling with larger maxDistance and '...
        'nUnitCell parameters!'],param.bond(1));
end

% select bonds with given atoms
% convert atom indices from the unit_cell into matom indices
if ~isempty(param.atom)
    
    matom = obj.matom;
    maIdx1a = ismember(obj.coupling.atom1,find(ismember(matom.idx,aIdx1)));
    maIdx2a = ismember(obj.coupling.atom2,find(ismember(matom.idx,aIdx2)));
    maIdx1b = ismember(obj.coupling.atom1,find(ismember(matom.idx,aIdx2)));
    maIdx2b = ismember(obj.coupling.atom2,find(ismember(matom.idx,aIdx1)));
    idx = idx & ((maIdx1a & maIdx2a) | (maIdx1b & maIdx2b));
end

idx = find(idx);

% check that all to be assigned bonds are symmetry generated
if obj.symmetry
    if obj.coupling.idx(max(idx)) > obj.coupling.nsym
        warning('spinw:addcoupling:SymmetryLimit','The assigned bonds are not generated using space group symmetry!');
    end
end

% if subIdx is defined, subselect bonds
if ~isempty(param.subIdx)
    idx = idx(param.subIdx);
end

if isempty(idx)
    error('spinw:addcoupling:NoBond','No matrix assigned, since no bond fulfilled the given conditions!')
end

Jmod   = obj.coupling.mat_idx(:,idx);
Tmod   = obj.coupling.type(:,idx);
Symmod = obj.coupling.sym(:,idx);

param.mat  = int32(param.mat);
param.type = int32(param.type);
param.sym  = int32(param.sym);

if any(ismember(Jmod(:),param.mat))
    warning('spinw:addcoupling:CouplingIdxWarning',['Same matrix already '...
        'assigned on some coupling, duplicate assigments are removed!']);
end

if isempty(param.mat)
    error('spinw:addcoupling:WrongMatrixLabel','The selected matrix label does not exists!')
end

if any(Jmod(3,:))
    error('spinw:addcoupling:TooManyCoupling',['The maximum '...
        'number of allowed couplings (3) per bond is reached!']);
end

for ii = 1:numel(param.mat)
    % generate one index for each column in Jmod
    idxSel = sub2ind(size(Jmod),sum(Jmod>0,1)+1,1:size(Jmod,2));
    % remove index where Jmod already contains param.mat
    idxSel = idxSel(~any(Jmod==param.mat(ii),1));
        
    Jmod(idxSel)   = param.mat(ii);
    Tmod(idxSel)   = param.type(ii);
    Symmod(idxSel) = param.sym(ii);
end

if obj.symmetry && nosympar && any(~param.sym)
    warning('spinw:addcoupling:SymetryLowered','By subselecting equivalent bonds, the symmetry of the corresponding bond(s) are reduced to P1!');
end

obj.coupling.mat_idx(:,idx) = Jmod;
obj.coupling.type(:,idx)    = Tmod;
obj.coupling.sym(:,idx)     = Symmod;

end