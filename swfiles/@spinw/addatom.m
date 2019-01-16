function addatom(obj, varargin)
% adds new atom
% 
% ### Syntax
% 
% `addatom(obj,Name,Value)`
% 
% ### Description
% 
% `addatom(obj,Name,Value)` adds a new atom to the list of symmetry
% inequivalent sites together with its properties, such as position, spin
% quantum number, form factor, etc.
% 
% ### Examples
% 
% To add a magnetic atom with $S=1$ at position $r=(0,0,0)$ and a
% non-magnetic one at $r=(1/2,0,0)$ with red and blue color respectively
% use the following command
%
% ```
% >>crystal = spinw;
% >>crystal.genlattice('lat_const',[4 3 3])
% >>crystal.addatom('r',[0 1/2; 0 0; 0 0],'S',[1 0],'color',{'red' 'blue'})
% >>crystal.plot
% ```
% 
% ### Input Arguments
% 
% `obj`
% : [spinw] object.
% 
% ### Name-Value Pair Arguments
% 
% `r`
% : Atomic positions stored in a matrix with dimensions of $[3\times
%   n_{atom}]$.
% 
% `label`
% : Names of the atoms in a cell for plotting and form factor
%   calculations (see [magion.dat]), e.g. `label={'atom1' 'atom2'
%   'atom3'}`.
%   Default value is `atomi`, where `i` is the atom index.
% 
% `S`
% : Spin quantum number stored in a row vector with $n_{atom}$ elements,
%   for non-magnetic atoms set S to zero. If not given the spin quantum
%   number is guessed from the given label of the atom. For example if
%   `label` is `MCr3+` or `Cr3+` then the $S=3/2$ high spin state is
%   assumed for Cr$^{3+}$. The spin values for every ion is stored in the
%   [magion.dat] file. If the atom type is unknown $S=0$ is assumed.
% 
% `color`
% : RGB color of the atoms for plotting stored in a matrix with dimensions
%   of $[3\times n_{atom}]$, where each column describes an RGB color. Each
%   value is between 0 and 255. Default value is the color stored in the
%   [atom.dat] file. Alternatively a name of the color can be given as a
%   string, for example `'White'`, for multiple atoms package it into a
%   cell. For the list of colors, see [swplot.color] or the [color.dat]
%   file.
% 
% `ox`
% : Oxidation number given as a double or it will be determined
%   automatically from label. Default value is 0.
% 
% `occ`
% : Occupancy, given as double. Default value is 1.
% 
% `formfact`
% : Neutron scattering form factor, given as a row vector with 9 numbers,
%   for details see [sw_mff]. Also string labels can be used from the
%   [magion.dat] file.
% 
% `formfactn`
% : Same as the `formfact` option.
% 
% `formfactx`
% : X-ray scattering form factor, given as 9 numbers, for details
%   see [sw_cff], also labels can be used from the [xrayion.dat] file.
% 
% `Z`
% : Atomic number, given as integer or determined from the atom label
%   automatically. Default value is 113 (Unobtanium).
% 
% `A`
% : Atomic mass, given as integer. Default is -1 for the natural
%   mixture of isotopes.
% 
% `bn`
% : Neutron scattering length, given as double. Not implemented yet.
% 
% `bx`
% : X-ray scattering length.
% 
% `biso`
% : Isotropic displacement factors in units of \\ang$^2$.
%   Definition is the same as in
%   [FullProf](https://www.ill.eu/sites/fullprof/), defining the
%   Debye-Waller factor as $W(d) = 1/8*b_{iso}/d^2$, which is included in
%   the structure factor as $exp(-2W(d))$.
% 
% `update`
% : If `true`, existing atom with the same label and position as a
%   new one will be updated. Default is `true`.
% 
% ### Output Arguments
% 
% The function modifies the [spinw.unit_cell] property of the obj
% [spinw] object.
% 
% ### See Also
% 
% [spinw.genlattice] \| [spinw.addmatrix] \| [swplot.color] \| [sw_mff] \| [sw_cff]
%

if nargin < 2
    swhelp spinw.addatom
    return
end

inpForm.fname  = {'r'     'S'      'label' 'color' 'ox'   'occ'  'bn'   'biso'  'update'};
inpForm.defval = {[]      []       {}      []      []     []     []     []      true    };
inpForm.size   = {[-1 -2] [-3 -4]  [1 -5] [1 -6]   [1 -7] [1 -7] [1 -7] [1 -10] [1 1]   };
inpForm.soft   = {false   true     true    true    true   true   true   true    false   };

inpForm.fname  = [inpForm.fname  {'bx'   'formfactn' 'formfactx' 'b'    'formfact' 'A'    'Z'   }];
inpForm.defval = [inpForm.defval {[]     []          []          []     []         []     []    }];
inpForm.size   = [inpForm.size   {[1 -7] [-1 -9]     [-1 -10]    [1 -7] [-8 -9]    [1 -7] [1 -7]}];
inpForm.soft   = [inpForm.soft   {true    true       true         true  true       true   true  }];

newAtom = sw_readparam(inpForm, varargin{:});

if isempty(newAtom.formfactn)
    newAtom.formfactn = newAtom.formfact;
end
if isempty(newAtom.bn)
    newAtom.bn = newAtom.b;
end

if ~any(size(newAtom.r)-[1 3])
    newAtom.r = newAtom.r';
end

% modulo
newAtom.r = mod(newAtom.r,1);

% save warnings
warn0 = warning;

if isempty(newAtom.label)
    % remove warnings since the user don't want specific atom
    warning('off','sw_mff:WrongInput');
    warning('off','sw_cff:WrongInput');
    warning('off','sw_nb:WrongInput');
%elseif ~any(newAtom.S)
elseif ~any(sw_sub1(newAtom.S))
    % the atom is not intentionally magnetic
    warning('off','sw_mff:WrongInput');
end

if size(newAtom.label,2) == 1
    newAtom.label = newAtom.label';
end

% number of old atoms
nOldAtom = size(obj.unit_cell.r,2);
% number of new atoms
nNewAtom = numel(newAtom.r)/3;

% isotropic displacement
if isempty(newAtom.biso)
    newAtom.biso = zeros(1,nNewAtom);
elseif numel(newAtom.biso)~= nNewAtom
    error('spinw:addatom:WrongInput','Wrong input for option ''biso''!');
end

% check oxidation number
if ~isempty(newAtom.ox)
    if numel(newAtom.ox) ~= nNewAtom
        error('spinw:addatom:WrongInput','Wrong input for option ''ox''!')
    end
end

% atomic number
if ~isempty(newAtom.Z)
    if numel(newAtom.Z) ~= nNewAtom
        error('spinw:addatom:WrongInput','Wrong input for option ''Z''!')
    end
end

% atomic mass number
if ~isempty(newAtom.A)
    if numel(newAtom.A) ~= nNewAtom
        error('spinw:addatom:WrongInput','Wrong input for option ''A''!')
    end
end

if  ~iscell(newAtom.label)
    newAtom.label = {newAtom.label};
end

if isempty(newAtom.occ)
    newAtom.occ = ones(1,nNewAtom);
end

if isempty(newAtom.A)
    newAtom.A = int32(-ones(1,nNewAtom));
end

if isempty(newAtom.label)
    % Generate atom labels if not given
    if ~isempty(newAtom.Z)
        newAtom.label = sw_atomdata(newAtom.Z,'name');
    else
        newAtom.label = repmat({'atom'},[1 nNewAtom]);
        newAtom.Z     = 113+zeros(1,nNewAtom); % Unobtanium
    end
    
    if isempty(newAtom.ox)
        newAtom.ox = zeros(1,nNewAtom);
    end
    for jj = 1:nNewAtom
        if newAtom.ox(jj)<0
            newAtom.label{jj} = [newAtom.label{jj} num2str(-newAtom.ox(jj)) '-_' num2str(jj+nOldAtom)];
        elseif newAtom.ox(jj)>0
            newAtom.label{jj} = [newAtom.label{jj} num2str(newAtom.ox(jj)) '+_' num2str(jj+nOldAtom)];
        else
            newAtom.label{jj} = [newAtom.label{jj} '_' num2str(jj+nOldAtom)];
        end
    end
elseif isempty(newAtom.Z) || isempty(newAtom.ox)
    % try to determine Z and ox values
    aOx  = zeros(1,nNewAtom);
    newAtom.name = cell(1,nNewAtom);
    
    for ii = 1:nNewAtom
        % remove the first word (good for typical .cif file labels)
        aName = strword(newAtom.label{ii},2,true);
        aName = aName{1};
        % get oxidation number
        aOx(ii) = str2double(aName(ismember(aName,'0':'9')));
        % sign
        if aName(end) == '-'
            aOx(ii) = -aOx(ii);
        end
        
        % remove numbers
        aName  = aName(ismember(aName,['A':'Z' 'a':'z']));
        % keep only the last 2 character
        if numel(aName)>2
            aName = aName(end+[-1 0]);
        end
        
        newAtom.name{ii} = aName;
    end
    
    aOx(isnan(aOx)) = 0;
    
    if isempty(newAtom.ox)
        newAtom.ox = aOx;
    end
    if isempty(newAtom.Z)
        newAtom.Z = sw_atomdata(newAtom.name,'Z');
    end
end

% Generate atom colors
if isempty(newAtom.color)
    newAtom.color = sw_atomdata(newAtom.Z,'color');
else
    newAtom.color = swplot.color(newAtom.color);
end

% Generate spins and magnetic form factor, default is generated from the label of the atom.
if size(newAtom.S,2) == 1
    newAtom.S = newAtom.S';
end

% get form factors

% get neutron scattering form factor if given
if isempty(newAtom.formfactn)
    [~,newAtom.ffn,newAtom.S0] = sw_mff(newAtom.label,[],11);
    % get the auto size of the magnetic moments if not given
    if isempty(newAtom.S)
        newAtom.S = newAtom.S0;
    end
end

if ischar(newAtom.formfactn)
    newAtom.formfactn = {newAtom.formfactn};
end

if isnumeric(newAtom.formfactn)
    newAtom.formfactn = newAtom.formfactn(:)';
    switch numel(newAtom.formfactn)
        case 1
            newAtom.formfactn = [zeros(1,10) newAtom.formfactn];
        case 9
            newAtom.formfactn = [newAtom.formfactn(1:8) 0 0 newAtom.formfactn(9)];
    end
end

if iscell(newAtom.formfactn)
    [~,newAtom.ffn] = sw_mff(newAtom.formfactn,[],11);
    %newAtom.ffn = permute(newAtom.ffn,[3 2 1]);
elseif ~isempty(newAtom.formfactn)
    newAtom.ffn = newAtom.formfactn;
end

newAtom.ffn = permute(newAtom.ffn,[3 2 1]);

% define non-magnetic atoms if S is empty
if isempty(newAtom.S)
    newAtom.S = zeros(1,nNewAtom);
end

% x-ray scattering form factor
if isempty(newAtom.formfactx)
    % 11 coefficients for the magnetic atoms
    [~,newAtom.ffx]            = sw_cff(newAtom.label);
end

if ischar(newAtom.formfactx)
    newAtom.formfactx = {newAtom.formfactx};
end
if iscell(newAtom.formfactx)
    [~,newAtom.ffx] = sw_cff(newAtom.formfactx);
    %newAtom.ffx = permute(newAtom.ffx,[3 2 1]);
elseif ~isempty(newAtom.formfactx)
    newAtom.ffx = newAtom.formfactx;
end

newAtom.ffx = permute(newAtom.ffx,[3 2 1]);

% include 2 zeros to make both form factor the same size
%newAtom.ffn = [newAtom.ffn(1,1:8,:) zeros(1,2,size(newAtom.ffn,3)) newAtom.ffn(1,9,:)];

newAtom.ff = [newAtom.ffn;newAtom.ffx];
newAtom.b  = ones(2,nNewAtom);

% get neutron scattering length
if isempty(newAtom.bn)
    newAtom.b(1,:) = sw_nb(newAtom.label);
else
    newAtom.b(1,:) = newAtom.bn;
end


newAtom.Z  = int32(newAtom.Z);


if obj.symbolic
    if ~isa(newAtom.S,'sym')
        symS = sym(0);
        for jj = 1:nNewAtom
            if newAtom.S(jj) > 0
                symS(jj) = sym(['S_' num2str(jj+nOldAtom)],'positive');
            else
                symS(jj) = sym(0);
            end
        end
        newAtom.S = symS;
    end
end

newAtom.color     = int32(newAtom.color);
newObj.unit_cell  = newAtom;
spinw.validate(newObj,'unit_cell');

cField = {'r' 'label' 'S' 'color' 'ox' 'occ' 'b' 'A' 'Z' 'biso'};

if newAtom.update
    % TODO check for structures with atoms having identical labels
    % find identical labels and positions between the new atoms and
    % existing atoms
    [iLabel, lIdx] = ismember(newObj.unit_cell.label,obj.unit_cell.label);
    [iPos, pIdx]   = ismember(newObj.unit_cell.r',obj.unit_cell.r','rows');
    % index in the new atoms
    newIdx  = iLabel&iPos';
    % index in the existing atoms
    oldIdx = lIdx(ismember(lIdx,pIdx));

    if any(newIdx)
        for ii = 3:numel(cField)
            % update values
            obj.unit_cell.(cField{ii})(:,oldIdx) = newObj.unit_cell.(cField{ii})(:,newIdx);
            % remove these new atoms
            newObj.unit_cell.(cField{ii}) = newObj.unit_cell.(cField{ii})(:,~newIdx);
        end
        % same for form factor
        obj.unit_cell.ff(:,:,oldIdx) = newObj.unit_cell.ff(:,:,newIdx);
        newObj.unit_cell.ff          = newObj.unit_cell.ff(:,:,~newIdx);
    end
end

if ~isempty(newObj.unit_cell.S)
    % if there are some unique new atoms left add them
    for ii = 1:numel(cField)
        obj.unit_cell.(cField{ii}) = [obj.unit_cell.(cField{ii}) newObj.unit_cell.(cField{ii})];
    end
    
    obj.unit_cell.ff    = cat(3,obj.unit_cell.ff,newObj.unit_cell.ff);
end
%spinw.validate(obj);

[~,~,rIdx] = unique(obj.unit_cell.r','rows');

% restore warnings
warning(warn0);

% check occupancy
if any(accumarray(rIdx,obj.unit_cell.occ)>1)
    warning('spinw:addatom:WrongInput','Occupancy on some site is larger than 1!')
end

end
