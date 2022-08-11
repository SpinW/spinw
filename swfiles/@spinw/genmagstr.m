function genmagstr(obj, varargin)
% generates magnetic structure
% 
% ### Syntax
% 
% `genmagstr(obj,Name,Value)`
% 
% ### Description
% 
% `genmagstr(obj,Name,Value)` is a Swiss knife for generating magnetic
% structures. It allows the definition of magnetic structures using several
% different ways, depending on the `mode` parameter. The generated magnetic
% structure is stored in the [obj.mag_str] field. The magetic structure is
% stored as Fourier components with arbitrary number of wave vectors in the
% [spinw] object. However spin waves can be only calculated if the magnetic
% structure has a single propagation vector (plus a k=0 ferromagnetic
% component perpendicular to the incommensurate magnetization), we simply
% call it single-k magnetic structure. Thus `genmagstr` enables to input
% magnetic structures that comply with this restriction by defining a
% magnetic structure by the moment directions (`S`) in the crystallographic
% cell, a propagation vector (`km`) and a vector that defines the normal of
% the rotation of the spin spiral (`n`). The function converts these values
% into Fourier components to store. To solve spin wave dispersion of
% magnetic structures with multiple propagation vectors, a magnetic
% supercell has to be defined where the propagation vector can be
% approximated to zero.
% 
% ### Examples
% 
% The example creates the multi-k magnetic structure of USb given by the
% `FQ` Fourier components and plots the magnetic structure:
% 
% ```
% >>USb = spinw
% >>USb.genlattice('lat_const',[6.203 6.203 6.203],'spgr','F m -3 m')
% >>USb.addatom('r',[0 0 0],'S',1)
% >>FQ = cat(3,[0;0;1+1i],[0;1+1i;0],[1+1i;0;0])>>
% >>k  = [0 0 1;0 1 0;1 0 0];
% >>USb.genmagstr('mode','fourier','S',FQ,'nExt',[1 1 1],'k',k)
% >>plot(USb,'range',[1 1 1])
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
% `'mode'`
% : Mode that determines how the magnetic structure is generated:
%   * `'random'` (optionally reads `k`, `n`, `nExt`)
%           generates a random structure in the structural cell if no
%           other arguments are specified here or previously in this spinw
%           object. If `nExt` is specified all spins in the supercell are
%           randomised. If `k` (and optionally `n`) is specified the
%           structure is then treated as incommensurate (similarly to the
%           `'helical'` mode but with randomized spins in the first cell).
%   * `'direct'` (reads `S`, optionally reads `k`, `nExt`)
%           direct input of the magnetic structure using the 
%           parameters of the single-k magnetic structure.
%   * `'tile'` (reads `S`, optionally reads `nExt`)
%           Simply extends the existing or input structure
%           (`S`) into a magnetic supercell by replicating it.
%           If no structure is stored in the [spinw] object a random
%           structure is generated automatically. If defined,
%           `S` is used as starting structure for extension
%           overwriting the stored structure. If the original
%           structure is already extended with other size, only the
%           moments in the crystallographic cell wil be replicated.
%           Magnetic ordering wavevector `k` will be set to zero. To
%           generate structure with non-zero k, use the `'helical'` or
%           `'direct'` option.
%   * `'helical'` (reads `S`, optionally reads `n`, `k`, `r0`, `nExt`, `epsilon`)
%           generates helical structure in a single cell or in a
%           supercell. In contrary to the `'extend'` option the
%           magnetic structure is not generated by replication but
%           by rotation of the moments using the following formula:
%
%     $\mathbf{S}^{gen}_i(\mathbf{r}) = R(2 \pi \mathbf{k_m} \cdot \mathbf{r})\cdot \mathbf{S}_i$
%
%     where $S_i$ has either a single moment or as many moments
%           as the number of magnetic atoms in the crystallographic
%           cell. In the first case $r$ denotes the atomic
%           positions, while for the second case $r$ denotes the
%           position of the origin of the cell.
%   * `'rotate'` (optionally reads `S`, `phi`, `phid`, `n`, `nExt`)
%           uniform rotation of all magnetic moments with a
%           `phi` angle around the given `n` vector. If
%           `phi=0`, all moments are rotated so, that the first
%           moment is parallel to `n` vector in case of
%           collinear structure or in case of planar structure
%           `n` defines the normal of the plane of the magnetic
%           moments.
%   * `'func'` (reads `func`, `x0`)
%           function that defines the parameters of the single-k
%           magnetic structure: moment vectors, propagation vector
%           and normal vector from arbitrary parameters in the
%           following form:
%     ```
%     [S, k, n] = @(x)func(S0, x)
%     ```  
%     where `S` is matrix with dimensions of $[3\times n_{magExt}]$. `k` is
%           the propagation vector in a 3-element row vector. `n` is the
%           normal vector of the spin rotation plane also 3-element row
%           vector. The default value for `func` is `@gm_spherical3d`. For planar
%           magnetic structure use `@gm_planar`. Only `func` and `x`
%           have to be defined for this mode.
%  * `'fourier'` (reads `S`, optionally reads `k`, `r0`, `nExt`, `epsilon`)
%           same as `'helical'`, except the `S` option is taken as the
%           Fourier components, thus if it contains real numbers, it will
%           generate a sinusoidally modulated structure instead of
%           a spiral.
% 
% `'phi'`
% : Angle of rotation of the magnetic moments in radian. Default
%   value is 0.
% 
% `'phid'`
% : Angle of rotation of the magnetic moments in \\deg. Default
%   value is 0.
% 
% `'nExt'`
% : Size of the magnetic supercell in multiples of the
%   crystallographic cell, dimensions are $[1\times 3]$. Default value is
%   stored in `obj`. If `nExt` is a single number, then the size of the
%   extended unit cell is automatically determined from the first
%   magnetic ordering wavevector. E.g. if `nExt = 0.01`, then the number
%   of unit cells is determined so, that in the extended unit cell,
%   the magnetic ordering wave vector is `[0 0 0]`, within the given
%   0.01 rlu tolerance.
% 
% `'k'`
% : Magnetic ordering wavevector in rlu, dimensions are $[n_K\times 3]$.
%   Default value is defined in `obj`.
% 
% `'n'`
% : Normal vector to the spin rotation plane for single-k magnetic
%   structures, stored in a 3-element row vector. Default value `[0 0 1]`. The
%   coordinate system of the vector is determined by `unit`.
% 
% `'S'`
% : Vector values of the spins (expectation value), dimensions are $[3\times n_{spin} n_K]$.
%   Every column defines the three $(S_x, S_y, S_z)$ components of
%   the spin (magnetic moment) in the $xyz$ Descartes coodinate system or
%   in lu. Default value is stored in `obj`.
% 
% `'unit'`
% : Determines the coordinate system for `S` and `n` vectors using a
%   string:
%   * `'xyz'`   Cartesian coordinate system, fixed to the lattice.
%               Default value.
%   * `'lu'`	Lattice coordinate system (not necessarily
%               Cartesian. The three coordinate vectors are
%               parallel to the lattice vectors but normalized to
%               unity.
% 
% `'epsilon'`
% : The smalles value of incommensurability that is
%   tolerated without warning in lattice units. Default is $10^{-5}$.
% 
% `'func'`
% : Function handle that produces the magnetic moments, ordering wave
%   vector and normal vector from the `x` parameters in the
%   following form:
%   ```
%   [M, k, n] = @(x)func(M0,x)
%   ```
%   where `M` is a matrix with dimensions of $[3\times n_{magExt}]$, `k` is
%   the propagation vector, `n` is the normal vector of the spin rotation
%   plane. The default function is [gm_spherical3d]. For planar magnetic
%   structure use [gm_planar].
% 
% `'x0'`
% : Input parameters for `func` function, row vector with $n_X$ number of
%   elements.
% 
% `'norm'`
% : Set the length of the generated magnetic moments to be equal to
%   the spin of the magnetic atoms. Default is `true`.
% 
% `'r0'`
% : If `true` and only a single spin direction is given, the spin
%   phases are determined by atomic position times k-vector, while
%   if it is `false`, the first spin will have 0 phase. Default is
%   `true`.
% 
% ### Output Arguments
% 
% The [obj.mag_str] field will contain the new magnetic structure.
% 
% ### See Also
% 
% [spinw] \| [spinw.anneal] \| [spinw.optmagstr] \| [gm_spherical3d] \| [gm_planar]
%
% *[rlu]: Reciprocal Lattice Units
% *[lu]: Lattice Units
%

if isempty(obj.matom.r)
    error('spinw:genmagstr:NoMagAtom','There are no magnetic atoms (S>0) in the unit cell!')
end

if isempty(obj.mag_str.k)
    k0 = [0 0 0];
else
    k0 = obj.mag_str.k';
end

inpForm.fname  = {'mode'   'nExt'            'k'           'n'     };
inpForm.defval = {'tile' obj.mag_str.nExt    []            nan(1,3)};
inpForm.size   = {[1 -1]   [1 -4]            [-6 3]        [-6 3]  };
inpForm.soft   = {false    false             true          false   };

inpForm.fname  = [inpForm.fname  {'func'          'x0'   'norm' 'r0' }];
inpForm.defval = [inpForm.defval {@gm_spherical3d []     true   true }];
inpForm.size   = [inpForm.size   {[1 1]           [1 -3] [1 1]  [1 1]}];
inpForm.soft   = [inpForm.soft   {false           true   false  false}];

inpForm.fname  = [inpForm.fname  {'S'       'phi' 'phid' 'epsilon' 'unit'}];
inpForm.defval = [inpForm.defval {[]         0     0      1e-5      'xyz'}];
inpForm.size   = [inpForm.size   {[3 -7 -6] [1 1] [1 1]  [1 1]     [1 -5]}];
inpForm.soft   = [inpForm.soft   {true      true  false  false     false }];

param = sw_readparam(inpForm, varargin{:});

% Error if S or k is provided but is empty. This means it has failed the
% input validation, but hasn't caused an error because soft=True
err_str = [];
if any(strcmp(varargin, 'S')) && isempty(param.S)
    err_str = ["S"];
end
if any(strcmp(varargin, 'k')) && isempty(param.k)
    err_str = [err_str "k"];
end
if length(err_str) > 0
    error('spinw:genmagstr:WrongInput', 'Incorrect input size for ' + ...
                                        join(err_str, ', '));
end

if strcmpi(param.mode, 'rotate') && isempty(obj.mag_str.F)
    error('spinw:genmagstr:WrongInput', ['rotate mode requires a magnetic ' ...
                                         'structure to be defined with another mode first'])
end

if isempty(param.k)
    noK = true;
    param.k = k0;
else
    noK = false;
end

if prod(double(param.nExt)) == 0
    error('spinw:genmagstr:WrongInput','''nExt'' has to be larger than 0!');
end

if strcmp(param.mode,'extend')
    param.mode = 'tile';
end

% input type for S, check whether it is complex type
cmplxS = ~isreal(param.S);
if strcmpi(param.mode, 'helical') && cmplxS
    error('spinw:genmagstr:WrongInput', ...
          ['S must be real for helical mode. To specify complex basis ' ...
           'vectors directly use fourier mode.'])
end

switch lower(param.unit)
    case 'lu'
        % convert the moments from lattice units to xyz
        BV = obj.basisvector(true);
        %param.S = BV*param.S;
        if ~isempty(param.S)
            param.S = mmat(BV,param.S);
        end
        param.n = (BV*param.n')';
    case 'xyz'
        % do nothing
    otherwise
        error('spinw:genmagstr:WrongInput','Option ''unit'' has to be either ''xyz'' or ''lu''!');
end

if isempty(param.S)
    % use the complex Fourier components from the stored magnetic structure
    param.S = obj.mag_str.F;
    cmplxS  = true;
end

% Magnetic ordering wavevector(s)
k  = param.k;
% number of k-vectors
nK = size(k,1);

nExt = double(param.nExt);

% automatic determination of the size of the extended unit cell based on
% the given k-vectors if nExt is a single number
if numel(nExt) == 1
    [~, nExtT] = rat(param.k(1,:),nExt);
    if nK>1
        for ii = 2:nK
            [~, nExtT2] = rat(param.k(ii,:),nExt);
            nExtT = lcm(nExtT,nExtT2);
        end
    end
    nExt = nExtT;
end

mAtom    = obj.matom;
nMagAtom = size(mAtom.r,2);
% number of magnetic atoms in the supercell
nMagExt  = nMagAtom*prod(nExt);

% Create mAtom.Sext matrix.
mAtom    = sw_extendlattice(nExt, mAtom);

% normalized axis of rotation, size (nK,3)
if isnan(param.n(1))
    % default value
    param.n = repmat([0 0 1],[nK 1]);
end
n = bsxfunsym(@rdivide,param.n,sqrt(sum(param.n.^2,2)));

if size(param.n,1) ~= nK
    error('spinw:genmagstr:WrongInput',['The number of normal vectors has'...
        ' to be equal to the number of k-vectors!'])
end

% convert input into symbolic variables
if obj.symb
    k       = sym(k);
    param.S = sym(param.S);
    n       = sym(n);
end

if ~cmplxS && ~strcmpi(param.mode,'fourier') && ~strcmpi(param.mode,'direct') && any(k(:))
    param.S = param.S + 1i*cross(repmat(permute(n,[2 3 1]),[1 size(param.S,2) 1]),param.S);
end

switch param.mode
    case 'tile'
        % effectively tiles the magnetic supercell with the given magnetic
        % moments if:
        % -the new number of extended cells does not equal to the number of
        %  cells defined in obj
        % -the number of spins stored in obj is not equal to the number
        %  of spins in the final structure
        if nMagAtom ~= size(param.S,2) && nMagExt ~= size(param.S,2)
            error('spinw:genmagstr:WrongInput', ['Incorrect input size for S, ' ...
                  'S must be provided for each magnetic atom']);
        elseif any(double(obj.mag_str.nExt) - double(param.nExt)) || (size(param.S,2) ~= nMagExt)
            S = repmat(param.S,[1 prod(nExt) 1]);
        else
            S = param.S;
        end
        % sum up all kvectors and keep the real part only
        S  = real(sum(S,3));
        k = [0 0 0];
        
    case 'random'
        % Create random spin directions and use a single k-vector
        S  = randn(nMagExt,3);
        S  = bsxfun(@rdivide,S,sqrt(sum(S.^2,2)));
        S  = bsxfunsym(@times,S,mAtom.Sext')';
        if noK
            k  = [0 0 0];
        end
        % keep the normal vector
        S = S + 1i*cross(repmat(permute(n,[2 3 1]),[1 size(S,2) 1]),S);
        
    case {'helical' 'fourier'}
        S0 = param.S;
        % Magnetic ordering wavevector in the extended unit cell.
        kExt = k.*nExt;
        % Warns about the non sufficient extension of the unit cell.
        % we substitute random values for symbolic km
        skExt = sw_sub1(kExt,'rand');
        if any(abs(skExt(:)-round(skExt(:)))>param.epsilon) && prod(nExt) > 1
            warning('spinw:genmagstr:UCExtNonSuff','In the extended unit cell k is still larger than epsilon!');
        end
        % number of spins in the input
        nSpin = size(param.S,2);
        
        if (nSpin~= nMagAtom) && (nSpin==1)
            % there is only a single given spin, use the fractional atomic position
            if param.r0
                r = mAtom.RRext;
            else
                r = bsxfun(@minus,mAtom.RRext,mAtom.RRext(:,1));
            end
            % add phase
            addPhase = true;
        elseif nSpin == nMagAtom
            % moments in the crystallographic unit cell are defined, use
            % only unit cell position.
            r = bsxfun(@rdivide,floor(bsxfun(@times,mAtom.RRext,nExt')),nExt');
            % add phase
            addPhase = true;

        elseif nSpin == nMagExt
            % no phase, magnetic structure is already in superstructure
            addPhase = false;
        else
            error('spinw:genmagstr:WrongNumberSpin','Wrong number of input spins!');
        end
               
        if addPhase
            % additional phase for each spin in the magnetic supercell
            %phi = sum(bsxfun(@times,2*pi*kExt',r),1);
            phi = sum(bsxfunsym(@times,2*pi*permute(kExt,[2 3 1]),r),1);
            
            % add the extra phase for each spin in the unit cell
            % TODO check
            S = bsxfunsym(@times,S0(:,mod(0:(nMagExt-1),nSpin)+1,:),exp(-1i*phi));
        else
            S = S0;
        end
    case 'direct'
        % direct input of real magnetic moments
        S = param.S;
        if size(S,2) == nMagAtom
            % single unit cell
            nExt = [1 1 1];
        end
        
        if size(S,2) ~= nMagExt
            error('spinw:genmagstr:WrongSpinSize','Wrong dimensions of param.S matrix!');
        end
                
    case 'rotate'
        
        if param.phi == 0
            % use degrees for phi if given
            param.phi = param.phid*pi/180;
        end
        
        S   = param.S;
        
        if ~isreal(param.phi)
            % rotate the first spin along [100]
            S1 = S(:,1)-sum(n*S(:,1))*n';
            S1 = S1/norm(S1);
            param.phi = -atan2(cross(n,[1 0 0])*S1,[1 0 0]*S1);
        end

        if param.phi == 0
            % The starting vector, size (1,3):
            incomm = mod(bsxfun(@times,k,nExt),1);
            incomm = any(incomm(:));
            if incomm
                S1 = sw_nvect(S);
            else
                S1 = sw_nvect(real(S));
            end
            
            % Axis of rotation defined by the spin direction
            nRot  = cross(n,S1);
            % Angle of rotation.
            phi = -atan2(norm(cross(S1,n)),dot(S1,n));
        else
            nRot = n;
            % Angle of rotation.
            phi = param.phi(1);
        end
        % Rotate the spins.
        S = sw_rot(nRot,phi,S);
        k = obj.mag_str.k';
        
    case 'func'
        S = mAtom.S;
        S = repmat(S,[prod(nExt) 1]);
        
        if obj.symbolic
            [S, k, n] = param.func(sym(S), sym(param.x0));
        else
            [S, k, n] = param.func(S,param.x0);
        end
        
        if any(k)
            S = S + 1i*cross(repmat(permute(n,[2 3 1]),[1 size(S,2) 1]),S);
        end
        
    otherwise
        error('spinw:genmagstr:WrongMode','Wrong param.mode value!');
end

% normalize the magnetic moments
if param.norm
    normS = sqrt(sum(real(S).^2,1))./repmat(mAtom.S,[1 prod(nExt)]);
    normS(normS==0) = 1;
    S = bsxfunsym(@rdivide,S,normS);
end

% simplify expressions
if obj.symbolic
    S = simplify(sym(S));
    k = simplify(sym(k));
end

mag_str.nExt = int32(nExt(:))';
mag_str.k    = k';
mag_str.F    = S;

obj.mag_str   = mag_str;
spinw.validate(obj);

end