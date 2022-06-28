function spectra = spinwave(obj, hkl, varargin)
% calculates spin correlation function using linear spin wave theory
%
% ### Syntax
%
% `spectra = spinwave(obj,Q)`
%
% `spectra = spinwave(___,Name,Value)`
%
% ### Description
%
% `spinwave(obj,Q,Name,Value)` calculates spin wave dispersion and
% spin-spin correlation function at the reciprocal space points $Q$. The
% function can solve any single-k magnetic structure exactly and any
% multi-k magnetic structure appoximately and quadratic spinw-spin
% interactions as well as single ion anisotropy and magnetic field.
% Biquadratic exchange interactions are also implemented, however only for
% $k_m=0$ magnetic structures.
%
% If the magnetic ordering wavevector is non-integer, the dispersion is
% calculated using a coordinate system rotating from unit cell to unit
% cell. In this case the spin Hamiltonian has to fulfill this extra
% rotational symmetry which is not checked programatically.
%
% Some of the code of the function can run faster if mex files are used. To
% switch on mex files, use the `swpref.setpref('usemex',true)` command. For
% details see the [sw_mex] and [swpref.setpref] functions.
%
% ### Examples
%
% To calculate and plot the spin wave dispersion of the
% triangular lattice antiferromagnet ($S=1$, $J=1$) along the $(h,h,0)$
% direction in reciprocal space we create the built in triangular lattice
% model using `sw_model`.
%
% ```
% >>tri = sw_model('triAF',1)
% >>spec = tri.spinwave({[0 0 0] [1 1 0]})
% >>sw_plotspec(spec)
% >>snapnow
% ```
% 
% ### Input Arguments
%
% `obj`
% : [spinw] object.
% 
% `Q`
% : Defines the $Q$ points where the spectra is calculated, in reciprocal
%   lattice units, size is $[3\times n_{Q}]$. $Q$ can be also defined by
%   several linear scan in reciprocal space. In this case `Q` is cell type,
%   where each element of the cell defines a point in $Q$ space. Linear scans
%   are assumed between consecutive points. Also the number of $Q$ points can
%   be specified as a last element, it is 100 by defaults. 
%   
%   For example to define a scan along $(h,0,0)$ from $h=0$ to $h=1$ using
%   200 $Q$ points the following input should be used:
%   ```
%   Q = {[0 0 0] [1 0 0]  200}
%   ```
%
%   For symbolic calculation at a general reciprocal space point use `sym`
%   type input. 
%
%   For example to calculate the spectrum along $(h,0,0)$ use:
%   ```
%   Q = [sym('h') 0 0]
%   ```
%   To calculate spectrum at a specific $Q$ point symbolically, e.g. at
%   $(0,1,0)$ use:
%   ```
%   Q = sym([0 1 0])
%   ```
%
% ### Name-Value Pair Arguments
%
% `'formfact'`
% : If true, the magnetic form factor is included in the spin-spin
%   correlation function calculation. The form factor coefficients are
%   stored in `obj.unit_cell.ff(1,:,atomIndex)`. Default value is `false`.
%
% `'formfactfun'`
% : Function that calculates the magnetic form factor for given $Q$ value.
%   value. Default value is `@sw_mff`, that uses a tabulated coefficients
%   for the form factor calculation. For anisotropic form factors a user
%   defined function can be written that has the following header:
%   ```
%   F = formfactfun(atomLabel,Q)
%   ```
%   where the parameters are:
%   * `F`           row vector containing the form factor for every input 
%                   $Q$ value
%   * `atomLabel`   string, label of the selected magnetic atom
%   * `Q`           matrix with dimensions of $[3\times n_Q]$, where each
%                   column contains a $Q$ vector in $\\ang^{-1}$ units.
%
% `'gtensor'`
% : If true, the g-tensor will be included in the spin-spin correlation
%   function. Including anisotropic g-tensor or different
%   g-tensor for different ions is only possible here. Including a simple
%   isotropic g-tensor is possible afterwards using the [sw_instrument]
%   function.
%
% `'fitmode'`
% : If `true`, function is optimized for multiple consecutive calls (e.g. 
%   the output spectrum won't contain the copy of `obj`), default is
%   `false`.
%
% `'notwin'`
% : If `true`, the spectra of the twins won't be calculated. Default is
% `false`.
%
% `'sortMode'`
% : If `true`, the spin wave modes will be sorted by continuity. Default is 
%   `true`.
%
% `'optmem'`
% : Parameter to optimise memory usage. The list of Q values will be cut
%   into `optmem` number of pieces and will be calculated piece by piece to
%   decrease peak memory usage. Default value is 0, when the number
%   of slices are determined automatically from the available free memory.
%
% `'tol'`
% : Tolerance of the incommensurability of the magnetic ordering wavevector.
%   Deviations from integer values of the ordering wavevector smaller than
%   the tolerance are considered to be commensurate. Default value is
%   $10^{-4}$.
%
% `'omega_tol'`
% : Tolerance on the energy difference of degenerate modes when
%   diagonalising the quadratic form, default value is $10^{-5}$.
%
% `'hermit'`
% : Method for matrix diagonalization with the following logical values:
% 
%   * `true`    using Colpa's method (for details see [J.H.P. Colpa, Physica 93A (1978) 327](http://www.sciencedirect.com/science/article/pii/0378437178901607)),
%               the dynamical matrix is converted into another Hermitian
%               matrix, that will give the real eigenvalues.
%   * `false`   using the standard method (for details see [R.M. White, PR 139 (1965) A450](https://journals.aps.org/pr/abstract/10.1103/PhysRev.139.A450))
%               the non-Hermitian $\mathcal{g}\times \mathcal{H}$ matrix
%               will be diagonalised, which is computationally less
%               efficient. Default value is `true`.
%
% {{note Always use Colpa's method, except when imaginary eigenvalues are
%   expected. In this case only White's method work. The solution in this
%   case is wrong, however by examining the eigenvalues it can give a hint
%   where the problem is.}}
%               
% `'saveH'`
% : If true, the quadratic form of the Hamiltonian is also saved in the
%   output. Be carefull, it can take up lots of memory. Default value is
%   `false`.
%
% `'saveV'`
% : If true, the matrices that transform the normal magnon modes into the
%   magnon modes localized on the spins are also saved into the output. Be
%   carefull, it can take up lots of memory. Default value is `false`.
%
% `'saveSabp'`
% : If true, the dynamical structure factor in the rotating frame
%   $S'(k,\omega)$ is saved. Default value is `false`.
%
% `'title'`
% : Gives a title string to the simulation that is saved in the output.
%
% `'fid'`
% : Defines whether to provide text output. The default value is determined
%   by the `fid` preference stored in [swpref]. The possible values are:
%   * `0`   No text output is generated.
%   * `1`   Text output in the MATLAB Command Window.
%   * `fid` File ID provided by the `fopen` command, the output is written
%           into the opened file stream.
%
% `'tid'`
% : Determines if the elapsed and required time for the calculation is
%   displayed. The default value is determined by the `tid` preference
%   stored in [swpref]. The following values are allowed (for more details
%   see [sw_timeit]):
%   * `0` No timing is executed.
%   * `1` Display the timing in the Command Window.
%   * `2` Show the timing in a separat pup-up window.
%
% ### Output Arguments
%
% `spectra`
% : structure, with the following fields:
%   * `omega`   Calculated spin wave dispersion with dimensions of
%               $[n_{mode}\times n_{Q}]$.
%   * `Sab`     Dynamical structure factor with dimensins of
%               $[3\times 3\times n_{mode}\times n_{Q}]$. Each
%               `(:,:,i,j)` submatrix contains the 9 correlation functions
%               $S^{xx}$, $S^{xy}$, $S^{xz}$, etc. If given, magnetic form
%               factor is included. Intensity is in \\hbar units, normalized
%               to the crystallographic unit cell.
%   * `H`       Quadratic form of the Hamiltonian. Only saved if `saveH` is
%               true.
%   * `V`       Transformation matrix from the normal magnon modes to the
%               magnons localized on spins using the following:
%               $x_i = \sum_j V_{ij} \times x_j'$
%               Only saved if `saveV` is true.
%   * `Sabp`    Dynamical structure factor in the rotating frame,
%               dimensions are $[3\times 3\times n_{mode}\times n_{Q}]$,
%               but the number of modes are equal to twice the number of
%               magnetic atoms.
%   * `formfact`  Cell containing the labels of the magnetic ions if form
%               factor in included in the spin-spin correlation function.
%   * `cmplxBase` The local coordinate system on each magnetic moment is
%               defined by the complex magnetic moments:
%               $\begin{align}  e_1 &= \Im(\hat{M})\\
%                               e_3 &= Re(\hat{M})\\
%                               e_2 &= e_3\times e_1
%               \end{align}$
%
%   * `hkl`     Contains the input $Q$ values, dimensions are $[3\times n_{Q}]$.
%   * `hklA`    Same $Q$ values, but in $\\ang^{-1}$ unit, in the
%               lab coordinate system, dimensins are $[3\times n_{Q}]$.
%   * `incomm`  Logical value, tells whether the calculated spectra is
%               incommensurate or not.
%   * `obj`     The copy (clone) of the input `obj`, see [spinw.copy].
%
% The number of magnetic modes (labeled by `nMode`) for commensurate
% structures is double the number of magnetic atoms in the magnetic cell.
% For incommensurate structures this number is tripled due to the
% appearance of the $(Q\pm k_m)$ Fourier components in the correlation
% functions. For every $Q$ points in the following order:
% $(Q-k_m,Q,Q+k_m)$.
%
% If several twins exist in the sample, `omega` and `Sab` are packaged into
% a cell, that contains $n_{twin}$ number of matrices.
%
% ### See Also
%
% [spinw] \| [spinw.spinwavesym] \| [sw_mex] \| [spinw.powspec] \| [sortmode]
%

pref = swpref;

% for linear scans create the Q line(s)
if nargin > 1
    hkl = sw_qscan(hkl);
else
    hkl = [];
end

% save warning of eigorth
orthWarn0 = false;

% save warning for singular matrix
singWarn0 = warning('off','MATLAB:nearlySingularMatrix');

% use mex file by default?
useMex = pref.usemex;

% calculate symbolic spectrum if obj is in symbolic mode
if obj.symbolic
    if numel(hkl) == 3
        hkl = sym(hkl);
    end
    
    if ~isa(hkl,'sym')
        inpForm.fname  = {'fitmode'};
        inpForm.defval = {false    };
        inpForm.size   = {[1 1]    };
        param0 = sw_readparam(inpForm, varargin{:});
        
        if ~param0.fitmode
            warning('spinw:spinwave:MissingInput','No hkl value was given, spin wave spectrum for general Q (h,k,l) will be calculated!');
        end
        spectra = obj.spinwavesym(varargin{:});
    else
        spectra = obj.spinwavesym(varargin{:},'hkl',hkl);
    end
    return
end

% help when executed without argument
if nargin==1
    swhelp spinw.spinwave
    spectra = [];
    return
end

title0 = 'Numerical LSWT spectrum';

inpForm.fname  = {'fitmode' 'notwin' 'sortMode' 'optmem' 'tol' 'hermit'};
inpForm.defval = {false     false    true       0        1e-4  true    };
inpForm.size   = {[1 1]     [1 1]    [1 1]      [1 1]    [1 1] [1 1]   };

inpForm.fname  = [inpForm.fname  {'omega_tol' 'saveSabp' 'saveV' 'saveH'}];
inpForm.defval = [inpForm.defval {1e-5        false      false   false  }];
inpForm.size   = [inpForm.size   {[1 1]       [1 1]      [1 1]   [1 1]  }];

inpForm.fname  = [inpForm.fname  {'formfact' 'formfactfun' 'title' 'gtensor'}];
inpForm.defval = [inpForm.defval {false       @sw_mff      title0  false    }];
inpForm.size   = [inpForm.size   {[1 -1]      [1 1]        [1 -2]  [1 1]    }];

inpForm.fname  = [inpForm.fname  {'cmplxBase' 'tid' 'fid' }];
inpForm.defval = [inpForm.defval {false       -1    -1    }];
inpForm.size   = [inpForm.size   {[1 1]       [1 1] [1 1] }];

param = sw_readparam(inpForm, varargin{:});

if ~param.fitmode
    % save the time of the beginning of the calculation
    spectra.datestart = datestr(now);
end

if param.fitmode
    param.sortMode = false;
    param.tid = 0;
end

if param.tid == -1
    param.tid = pref.tid;
end

if param.fid == -1
    param.fid = pref.fid;
end
fid = param.fid;

% generate magnetic structure in the rotating noation
magStr = obj.magstr;

% size of the extended magnetic unit cell
nExt    = magStr.N_ext;
% magnetic ordering wavevector in the extended magnetic unit cell
km = magStr.k.*nExt;

% whether the structure is incommensurate
incomm = any(abs(km-round(km)) > param.tol);

% Transform the momentum values to the new lattice coordinate system
hkl = obj.unit.qmat*hkl;

% Calculates momentum transfer in A^-1 units.
hklA = 2*pi*(hkl'/obj.basisvector)';

% Check for 2*km
tol = param.tol*2;
helical =  sum(abs(mod(abs(2*km)+tol,1)-tol).^2) > tol;

% number of Q points
nHkl0 = size(hkl,2);

% define Q scans for the twins
nTwin = size(obj.twin.vol,2);
if param.notwin
    nTwin = 1;
end

% if the single twin has no rotation set param.notwin true
rotc1 = obj.twin.rotc(:,:,1)-eye(3);
if (nTwin == 1) && norm(rotc1(:))==0
    param.notwin = true;
end

if ~param.notwin
    % In the abc coordinate system of the selected twin the scan is
    % rotated opposite direction to rotC.
    hkl  = obj.twinq(hkl);
    nHkl = nHkl0*nTwin;
else
    nHkl = nHkl0;
    hkl  = {hkl};
end

if incomm
    % TODO
    if ~helical && ~param.fitmode
        warning('spinw:spinwave:Twokm',['The two times the magnetic ordering '...
            'wavevector 2*km = G, reciproc lattice vector, use magnetic supercell to calculate spectrum!']);
    end
    
    hkl0 = cell(1,nTwin);
    hklExt = cell(1,nTwin);
    
    for tt = 1:nTwin
        % without the k_m: (k, k, k)
        hkl0{tt} = repmat(hkl{tt},[1 3]);
        
        % for wavevectors in the extended unit cell km won't be multiplied by
        % nExt (we devide here to cancel the multiplication later)
        kme = km./nExt;
        hklExt{tt}  = [bsxfun(@minus,hkl{tt},kme') hkl{tt} bsxfun(@plus,hkl{tt},kme')];
        
        % calculate dispersion for (k-km, k, k+km)
        hkl{tt}  = [bsxfun(@minus,hkl{tt},km') hkl{tt} bsxfun(@plus,hkl{tt},km')];
    end
    nHkl  = nHkl*3;
    nHkl0 = nHkl0*3;
else
    hkl0   = hkl;
    hklExt = hkl;
end

hkl    = cell2mat(hkl);
hkl0   = cell2mat(hkl0);
hklExt = cell2mat(hklExt);

% determines a twin index for every q point
twinIdx = repmat(1:nTwin,[nHkl0 1]);
twinIdx = twinIdx(:);

% Create the interaction matrix and atomic positions in the extended
% magnetic unit cell.
[SS, SI, RR] = obj.intmatrix('fitmode',true,'conjugate',true);

% add the dipolar interactions to SS.all
SS.all = [SS.all SS.dip];

% is there any biquadratic exchange
bq = SS.all(15,:)==1;

% Biquadratic exchange only supported for commensurate structures
if incomm && any(bq)
    error('spinw:spinwave:Biquadratic','Biquadratic exchange can be only calculated for k=0 structures!');
end

if any(bq)
    % Separate the biquadratic couplings
    % Just use the SS.bq matrix produced by intmatrix(), it won't contain
    % the transpose matrices (not necessary for biquadratic exchange)
    % TODO check whether to keep the transposed matrices to be sure
    SS.bq = SS.all(1:6,bq);
    % Keep only the quadratic exchange couplings
    SS.all = SS.all(1:14,SS.all(15,:)==0);
end

% Converts wavevctor list into the extended unit cell
hklExt  = bsxfun(@times,hklExt,nExt')*2*pi;
% q values without the +/-k_m value
hklExt0 = bsxfun(@times,hkl0,nExt')*2*pi;

% Calculates parameters eta and zed.
if isempty(magStr.S)
    error('spinw:spinwave:NoMagneticStr','No magnetic structure defined in obj!');
end

M0 = magStr.S;
S0 = sqrt(sum(M0.^2,1));
% normal to rotation of the magnetic moments
n  = magStr.n;
nMagExt = size(M0,2);

if incomm
    fprintf0(fid,['Calculating INCOMMENSURATE spin wave spectra '...
        '(nMagExt = %d, nHkl = %d, nTwin = %d)...\n'],nMagExt, nHkl0, nTwin);
else
    fprintf0(fid,['Calculating COMMENSURATE spin wave spectra '...
        '(nMagExt = %d, nHkl = %d, nTwin = %d)...\n'],nMagExt, nHkl0, nTwin);
end

% If cmplxBase is false, we use a local (e1,e2,e3) coordinate system fixed
% to the moments:
% e3||Si,ata
% e2 = Si x [1,0,0], if Si || [1,0,0] --> e2 = [0,0,1]
% e1 = e2 x e3
% If cmplxBase is true, we use a coordinate system fixed by the
% complex magnetisation vectors:
% e1 = imag(M)
% e3 = real(M)
% e2 = e3 x e1
if ~param.cmplxBase
    if obj.symbolic
        e3 = simplify(M0./[S0; S0; S0]);
        % e2 = Si x [1,0,0], if Si || [1,0,0] --> e2 = [0,0,1]
        e2  = [zeros(1,nMagExt); e3(3,:); -e3(2,:)];
        % select zero vector and make them parallel to [0,0,1]
        selidx = abs(e2)>0;
        if isa(selidx,'sym')
            e2(3,~any(~sw_always(abs(e2)==0))) = 1;
        else
            e2(3,~any(abs(e2)>0)) = 1;
        end
        E0 = sqrt(sum(e2.^2,1));
        e2  = simplify(e2./[E0; E0; E0]);
        % e1 = e2 x e3
        e1  = simplify(cross(e2,e3));
    else
        % e3 || Si
        e3 = bsxfun(@rdivide,M0,S0);
        % e2 = Si x [1,0,0], if Si || [1,0,0] --> e2 = [0,0,1]
        e2  = [zeros(1,nMagExt); e3(3,:); -e3(2,:)];
        e2(3,~any(abs(e2)>1e-10)) = 1;
        e2  = bsxfun(@rdivide,e2,sqrt(sum(e2.^2,1)));
        % e1 = e2 x e3
        e1  = cross(e2,e3);
    end
else
    F0  = obj.mag_str.F;
    RF0 = sqrt(sum(real(F0).^2,1));
    IF0 = sqrt(sum(imag(F0).^2,1));
    % e3 = real(M)
    e3  = real(F0)./repmat(RF0,[3 1]);
    % e1 = imag(M) perpendicular to e3
    e1  = imag(F0)./repmat(IF0,[3 1]);
    e1  = e1-bsxfun(@times,sum(e1.*e3,1),e3);
    e1  = e1./repmat(sqrt(sum(e1.^2,1)),[3 1]);
    % e2 = cross(e3,e1)
    e2  = cross(e3,e1);
    
    if obj.symbolic
        e1 = simplify(e1);
        e2 = simplify(e2);
        e3 = simplify(e3);
    end
    
end
% assign complex vectors that define the rotating coordinate system on
% every magnetic atom
zed = e1 + 1i*e2;
eta = e3;

dR    = [SS.all(1:3,:) zeros(3,nMagExt)];
atom1 = [SS.all(4,:)   1:nMagExt];
atom2 = [SS.all(5,:)   1:nMagExt];
% magnetic couplings, 3x3xnJ
JJ = cat(3,reshape(SS.all(6:14,:),3,3,[]),SI.aniso);

if incomm
    % transform JJ due to the incommensurate wavevector
    [~, K] = sw_rot(n,km*dR*2*pi);
    % multiply JJ with K matrices for every interaction
    % and symmetrising JJ for the rotating basis
    JJ = (mmat(JJ,K)+mmat(K,JJ))/2;
end

nCoupling = size(JJ,3);

zedL = repmat(permute(zed(:,atom1),[1 3 2]),[1 3 1]);
zedR = repmat(permute(zed(:,atom2),[3 1 2]),[3 1 1]);

etaL = repmat(permute(eta(:,atom1),[1 3 2]),[1 3 1]);
etaR = repmat(permute(eta(:,atom2),[3 1 2]),[3 1 1]);

SiSj = sqrt(S0(atom1).*S0(atom2));

% Creates temporary values for calculating matrix elements.
AD  =  shiftdim(sum(sum(etaL.*JJ.*etaR,2),1),1);
A20 = -S0(atom2).*AD;
D20 = -S0(atom1).*AD;
BC0 =  SiSj.*shiftdim(sum(sum(zedL.*JJ.*     zedR ,2),1),1);
AD0 =  SiSj.*shiftdim(sum(sum(zedL.*JJ.*conj(zedR),2),1),1);

% Magnetic field is different for every twin
%MF  =  repmat(obj.unit.muB*SI.field*eta,[1 2]);
MF = zeros(1,2*nMagExt,nTwin);
for ii = 1:nTwin
    % rotate the magnetic field to the relative direction of every twin
    % backward rotation with the rotc matrix of the twin
    twinB = SI.field*obj.twin.rotc(:,:,ii)*obj.unit.muB;
    MF(:,:,ii) = repmat(twinB*permute(mmat(SI.g,permute(eta,[1 3 2])),[1 3 2]),[1 2]);
end

% Creates the serial indices for every matrix element in ham matrix.
idxA1 = [atom1'         atom2'         ];
idxA2 = [atom1'         atom1'         ];
idxB  = [atom1'         atom2'+nMagExt ];
% transpose of idxB
%idxC  = [atom2'+nMagExt atom1'         ]; % SP1
idxD1 = idxA1+nMagExt;
idxD2 = [atom2'+nMagExt atom2'+nMagExt ];
idxMF = [(1:2*nMagExt)' (1:2*nMagExt)' ];

% Calculate matrix elements for biquadratic exchange
if any(bq)
    bqdR    = SS.bq(1:3,:);
    bqAtom1 = SS.bq(4,:);
    bqAtom2 = SS.bq(5,:);
    bqJJ    = SS.bq(6,:);
    nbqCoupling = numel(bqJJ);
    
    % matrix elements: M,N,P,Q
    bqM = sum(eta(:,bqAtom1).*eta(:,bqAtom2),1);
    bqN = sum(eta(:,bqAtom1).*zed(:,bqAtom2),1);
    bqO = sum(zed(:,bqAtom1).*zed(:,bqAtom2),1);
    bqP = sum(conj(zed(:,bqAtom1)).*zed(:,bqAtom2),1);
    bqQ = sum(zed(:,bqAtom1).*eta(:,bqAtom2),1);
    
    Si = S0(bqAtom1);
    Sj = S0(bqAtom2);
    % C_ij matrix elements
    bqA0 = (Si.*Sj).^(3/2).*(bqM.*conj(bqP) + bqQ.*conj(bqN)).*bqJJ;
    bqB0 = (Si.*Sj).^(3/2).*(bqM.*bqO + bqQ.*bqN).*bqJJ;
    bqC  = Si.*Sj.^2.*(conj(bqQ).*bqQ - 2*bqM.^2).*bqJJ;
    bqD  = Si.*Sj.^2.*(bqQ).^2.*bqJJ;
    
    % Creates the serial indices for every matrix element in ham matrix.
    % Aij(k) matrix elements (b^+ b)
    idxbqA  = [bqAtom1' bqAtom2'];
    % b b^+ elements
    idxbqA2 = [bqAtom1' bqAtom2']+nMagExt;
    
    % Bij(k) matrix elements (b^+ b^+)
    idxbqB  = [bqAtom1' bqAtom2'+nMagExt];
    % transpose of B (b b)
    %idxbqB2 = [bqAtom2'+nMagExt bqAtom1']; % SP2
    
    idxbqC  = [bqAtom1' bqAtom1'];
    idxbqC2 = [bqAtom1' bqAtom1']+nMagExt;
    
    idxbqD  = [bqAtom1' bqAtom1'+nMagExt];
    %idxbqD2 = [bqAtom1'+nMagExt bqAtom1]; % SP2
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEMORY MANAGEMENT LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if param.optmem == 0
    freeMem = sw_freemem;
    if freeMem > 0
        nSlice = ceil(nMagExt^2*nHkl*6912/freeMem*2);
    else
        nSlice = 1;
        if ~param.fitmode
            warning('spinw:spinwave:FreeMemSize','The size of the free memory is unkown, no memory optimisation!');
        end
    end
else
    nSlice = param.optmem;
end

if nHkl < nSlice
    fprintf0(fid,['Memory allocation is not optimal, nMagExt is'...
        ' too large compared to the free memory!\n']);
    nSlice = nHkl;
elseif nSlice > 1
    fprintf0(fid,['To optimise memory allocation, Q is cut'...
        ' into %d pieces!\n'],nSlice);
end

% message for magnetic form factor calculation
yesNo = {'No' 'The'};
fprintf0(fid,[yesNo{param.formfact+1} ' magnetic form factor is'...
    ' included in the calculated structure factor.\n']);
% message for g-tensor calculation
fprintf0(fid,[yesNo{param.gtensor+1} ' g-tensor is included in the '...
    'calculated structure factor.\n']);

if param.gtensor
    
    gtensor = SI.g;
    
    if incomm
        % keep the rotation invariant part of g-tensor
        nx  = [0 -n(3) n(2);n(3) 0 -n(1);-n(2) n(1) 0];
        nxn = n'*n;
        m1  = eye(3);
        gtensor = 1/2*gtensor - 1/2*mmat(mmat(nx,gtensor),nx) + 1/2*mmat(mmat(nxn-m1,gtensor),nxn) + 1/2*mmat(mmat(nxn,gtensor),2*nxn-m1);
    end
end

hklIdx = [floor(((1:nSlice)-1)/nSlice*nHkl)+1 nHkl+1];

% Empty omega dispersion of all spin wave modes, size: 2*nMagExt x nHkl.
omega = zeros(2*nMagExt, nHkl);

% empty Sab
Sab = zeros(3,3,2*nMagExt,0);

% Empty matrices to save different intermediate results for further
% analysis: Hamiltonian, eigenvectors, dynamical structure factor in the
% rotating frame
if param.saveV
    Vsave = zeros(2*nMagExt,2*nMagExt,nHkl);
end
if param.saveH
    Hsave = zeros(2*nMagExt,2*nMagExt,nHkl);
end

sw_timeit(0,1,param.tid,'Spin wave spectrum calculation');

warn1 = false;

% calculate all magnetic form factors
if param.formfact
    spectra.formfact = true;
    % Angstrom^-1 units for Q
    hklA0 = 2*pi*(hkl0'/obj.basisvector)';
    % store form factor per Q point for each atom in the magnetic supercell
    % TODO check prod(nExt)? instead of nExt
    %FF = repmat(param.formfactfun(permute(obj.unit_cell.ff(1,:,obj.matom.idx),[3 2 1]),hklA0),[1 nExt]);
    FF = repmat(param.formfactfun(permute(obj.unit_cell.ff(1,:,obj.matom.idx),[3 2 1]),hklA0),[prod(nExt) 1]);
else
    spectra.formfact = false;
end

for jj = 1:nSlice
    % q indices selected for every chunk
    hklIdxMEM  = hklIdx(jj):(hklIdx(jj+1)-1);
    % q values contatining the k_m vector
    hklExtMEM  = hklExt(:,hklIdxMEM);
    % q values without the +/-k_m vector
    hklExt0MEM = hklExt0(:,hklIdxMEM);
    % twin indices for every q point
    twinIdxMEM = twinIdx(hklIdxMEM);
    nHklMEM = size(hklExtMEM,2);
        
    % Creates the matrix of exponential factors nCoupling x nHkl size.
    % Extends dR into 3 x 3 x nCoupling x nHkl
    %     ExpF = exp(1i*permute(sum(repmat(dR,[1 1 nHklMEM]).*repmat(...
    %         permute(hklExtMEM,[1 3 2]),[1 nCoupling 1]),1),[2 3 1]))';
    ExpF = exp(1i*permute(sum(bsxfun(@times,dR,permute(hklExtMEM,[1 3 2])),1),[2 3 1]))';
    
    % Creates the matrix elements containing zed.
    A1 = bsxfun(@times,     AD0 ,ExpF);
    B  = bsxfun(@times,     BC0 ,ExpF);
    D1 = bsxfun(@times,conj(AD0),ExpF);
    
    
    
    % Store all indices
    % SP1: speedup for creating the matrix elements
    %idxAll = [idxA1; idxB; idxC; idxD1]; % SP1
    idxAll   = [idxA1; idxB; idxD1];
    % Store all matrix elements
    %ABCD   = [A1     B     conj(B)  D1]; % SP1
    ABCD   = [A1     2*B      D1];
    
    % Stores the matrix elements in ham.
    %idx3   = repmat(1:nHklMEM,[4*nCoupling 1]); % SP1
    idx3   = repmat(1:nHklMEM,[3*nCoupling 1]);
    idxAll = [repmat(idxAll,[nHklMEM 1]) idx3(:)];
    idxAll = idxAll(:,[2 1 3]);
    
    ABCD   = ABCD';
    
    
    % quadratic form of the boson Hamiltonian stored as a square matrix
    ham = accumarray(idxAll,ABCD(:),[2*nMagExt 2*nMagExt nHklMEM]);
    
    ham = ham + repmat(accumarray([idxA2; idxD2],2*[A20 D20],[1 1]*2*nMagExt),[1 1 nHklMEM]);
    
    if any(bq)
        % bqExpF = exp(1i*permute(sum(repmat(bqdR,[1 1 nHklMEM]).*repmat(...
        %     permute(hklExtMEM,[1 3 2]),[1 nbqCoupling 1]),1),[2 3 1]))';
        bqExpF = exp(1i*permute(sum(bsxfun(@times,bqdR,permute(hklExtMEM,[1 3 2])),1),[2 3 1]))';
        
        bqA  = bsxfun(@times,     bqA0, bqExpF);
        bqA2 = bsxfun(@times,conj(bqA0),bqExpF);
        bqB  = bsxfun(@times,     bqB0, bqExpF);
        idxbqAll = [idxbqA; idxbqA2; idxbqB];
        %bqABCD = [bqA bqA2 2*bqB];
        bqABCD = [bqA bqA2 2*bqB];
        bqidx3   = repmat(1:nHklMEM,[3*nbqCoupling 1]);
        idxbqAll = [repmat(idxbqAll,[nHklMEM 1]) bqidx3(:)];
        idxbqAll = idxbqAll(:,[2 1 3]);
        bqABCD = bqABCD';
        % add biquadratic exchange
        ham = ham + accumarray(idxbqAll,bqABCD(:),[2*nMagExt 2*nMagExt nHklMEM]);
        % add diagonal terms
        ham = ham + repmat(accumarray([idxbqC; idxbqC2; idxbqD],[bqC bqC 2*bqD],[1 1]*2*nMagExt),[1 1 nHklMEM]);
        
    end
    if any(SI.field)
        % different field for different twin
        for ii = min(twinIdxMEM):max(twinIdxMEM)
            nTwinQ = sum(twinIdxMEM==ii);
            ham(:,:,twinIdxMEM==ii) = ham(:,:,twinIdxMEM==ii) + ...
                repmat(accumarray(idxMF,MF(:,:,ii),[1 1]*2*nMagExt),[1 1 nTwinQ]);
        end
        
        %ham = ham + repmat(accumarray(idxMF,MF,[1 1]*2*nMagExt),[1 1 nHklMEM]);
    end
    
    ham = (ham + conj(permute(ham,[2 1 3])))/2;
    
    % diagonal of the boson commutator matrix
    gCommd = [ones(nMagExt,1); -ones(nMagExt,1)];
    % boson commutator matrix
    gComm  = diag(gCommd);
    %gd = diag(g);
    
    if param.hermit
        % All the matrix calculations are according to Colpa's paper
        % J.H.P. Colpa, Physica 93A (1978) 327-353
        
        % basis functions of the magnon modes
        V = zeros(2*nMagExt,2*nMagExt,nHklMEM);
        
        if useMex && nHklMEM>1
            % use mex files to speed up the calculation
            % mex file will return an error if the matrix is not positive definite.
            [K2, invK] = chol_omp(ham,'Colpa','tol',param.omega_tol);
            [V, omega(:,hklIdxMEM)] = eig_omp(K2,'sort','descend');
            % the inverse of the para-unitary transformation V
            for ii = 1:nHklMEM
                V(:,:,ii) = V(:,:,ii)*diag(sqrt(gCommd.*omega(:,hklIdxMEM(ii))));
            end
            % V = bsxfun(@times,invK,V);
            V = sw_mtimesx(invK,V);
        else
            for ii = 1:nHklMEM
                [K, posDef]  = chol(ham(:,:,ii));
                if posDef > 0
                    try
                        % get tolerance from smallest negative eigenvalue
                        tol0 = eig(ham(:,:,ii));
                        tol0 = sort(real(tol0));
                        tol0 = abs(tol0(1));
                        % TODO determine the right tolerance value
                        tol0 = tol0*sqrt(nMagExt*2)*4;
                        if tol0>param.omega_tol
                            error('spinw:spinwave:NonPosDefHamiltonian','Very baaaad!');
                        end
                        try
                            K = chol(ham(:,:,ii)+eye(2*nMagExt)*tol0);
                        catch
                            K = chol(ham(:,:,ii)+eye(2*nMagExt)*param.omega_tol);
                        end
                        warn1 = true;
                    catch PD
                        if param.tid == 2
                            % close timer window
                            sw_timeit(100,2,param.tid);
                        end
                        error('spinw:spinwave:NonPosDefHamiltonian',...
                            ['Hamiltonian matrix is not positive definite, probably'...
                            ' the magnetic structure is wrong! For approximate'...
                            ' diagonalization try the param.hermit=false option']);
                    end
                end
                
                K2 = K*gComm*K';
                K2 = 1/2*(K2+K2');
                % Hermitian K2 will give orthogonal eigenvectors
                [U, D] = eig(K2);
                D      = diag(D);
                
                % sort modes accordign to the real part of the energy
                [~, idx] = sort(real(D),'descend');
                U = U(:,idx);
                % omega dispersion
                omega(:, hklIdxMEM(ii)) = D(idx);
                
                % the inverse of the para-unitary transformation V
                V(:,:,ii) = inv(K)*U*diag(sqrt(gCommd.*omega(:, hklIdxMEM(ii)))); %#ok<MINV>
            end
        end
    else
        % All the matrix calculations are according to White's paper
        % R.M. White, et al., Physical Review 139, A450?A454 (1965)
        
        gham = mmat(gComm,ham);
        
        [V, D, orthWarn] = eigorth(gham,param.omega_tol,useMex);
        
        orthWarn0 = orthWarn || orthWarn0;
        
        for ii = 1:nHklMEM
            % multiplication with g removed to get negative and positive
            % energies as well
            omega(:,hklIdxMEM(ii)) = D(:,ii);
            M              = diag(gComm*V(:,:,ii)'*gComm*V(:,:,ii));
            V(:,:,ii)      = V(:,:,ii)*diag(sqrt(1./M));
        end
    end
    
    if param.saveV
        Vsave(:,:,hklIdxMEM) = V;
    end
    if param.saveH
        Hsave(:,:,hklIdxMEM) = ham;
    end
    
    % Calculates correlation functions.
    % V right
    VExtR = repmat(permute(V  ,[4 5 1 2 3]),[3 3 1 1 1]);
    % V left: conjugate transpose of V
    VExtL = conj(permute(VExtR,[1 2 4 3 5]));
    
    % Introduces the exp(-ikR) exponential factor.
    ExpF =  exp(-1i*sum(repmat(permute(hklExt0MEM,[1 3 2]),[1 nMagExt 1]).*repmat(RR,[1 1 nHklMEM]),1));
    % Includes the sqrt(Si/2) prefactor.
    ExpF = ExpF.*repmat(sqrt(S0/2),[1 1 nHklMEM]);
    
    ExpFL =      repmat(permute(ExpF,[1 4 5 2 3]),[3 3 2*nMagExt 2]);
    % conj transpose of ExpFL
    ExpFR = conj(permute(ExpFL,[1 2 4 3 5]));
    
    zeda = repmat(permute([zed conj(zed)],[1 3 4 2]),[1 3 2*nMagExt 1 nHklMEM]);
    % conj transpose of zeda
    zedb = conj(permute(zeda,[2 1 4 3 5]));
    
    % calculate magnetic structure factor using the hklExt0 Q-values
    % since the S(Q+/-k,omega) correlation functions also belong to the
    % F(Q)^2 form factor
    
    if param.formfact
        % include the form factor in the z^alpha, z^beta matrices
        zeda = zeda.*repmat(permute(FF(:,hklIdxMEM),[3 4 5 1 2]),[3 3 2*nMagExt 2 1]);
        zedb = zedb.*repmat(permute(FF(:,hklIdxMEM),[3 4 1 5 2]),[3 3 2 2*nMagExt 1]);
    end
    
    if param.gtensor
        % include the g-tensor
        zeda = mmat(repmat(permute(gtensor,[1 2 4 3]),[1 1 1 2]),zeda);
        zedb = mmat(zedb,repmat(gtensor,[1 1 2]));
    end
    % Dynamical structure factor from S^alpha^beta(k) correlation function.
    % Sab(alpha,beta,iMode,iHkl), size: 3 x 3 x 2*nMagExt x nHkl.
    % Normalizes the intensity to single unit cell.
    Sab = cat(4,Sab,squeeze(sum(zeda.*ExpFL.*VExtL,4)).*squeeze(sum(zedb.*ExpFR.*VExtR,3))/prod(nExt));
    
    sw_timeit(jj/nSlice*100,0,param.tid);
end

[~,singWarn] = lastwarn;
% restore warning for singular matrix
warning(singWarn0.state,'MATLAB:nearlySingularMatrix');

% If number of formula units are given per cell normalize to formula
% unit
if obj.unit.nformula > 0
    Sab = Sab/double(obj.unit.nformula);
end

sw_timeit(100,2,param.tid);

fprintf0(fid,'Calculation finished.\n');

if warn1 && ~param.fitmode
    warning('spinw:spinwave:NonPosDefHamiltonian',['To make the Hamiltonian '...
        'positive definite, a small omega_tol value was added to its diagonal!'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END MEMORY MANAGEMENT LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if incomm
    % resize matrices due to the incommensurability (k-km,k,k+km) multiplicity
    kmIdx = repmat(sort(repmat([1 2 3],1,nHkl0/3)),1,nTwin);
    % Rodrigues' rotation formula.
    nx  = [0 -n(3) n(2); n(3) 0 -n(1); -n(2) n(1) 0];
    nxn = n'*n;
    K1 = 1/2*(eye(3) - nxn - 1i*nx);
    K2 = nxn;
    
    % keep the rotation invariant part of Sab
    %nx  = [0 -n(3) n(2);n(3) 0 -n(1);-n(2) n(1) 0];
    %nxn = n'*n;
    m1  = eye(3);
    
    % if the 2*km vector is integer, the magnetic structure is not a true
    % helix
    %tol = param.tol*2;
    %helical =  sum(abs(mod(abs(2*km)+tol,1)-tol).^2) > tol;
    
    if helical
        % integrating out the arbitrary initial phase of the helix
        Sab = 1/2*Sab - 1/2*mmat(mmat(nx,Sab),nx) + 1/2*mmat(mmat(nxn-m1,Sab),nxn) + 1/2*mmat(mmat(nxn,Sab),2*nxn-m1);
    end
    
    % Save the structure factor in the rotating frame
    if param.saveSabp
        Sabp = Sab(:,:,:,kmIdx==2);
        omegap = omega(:,kmIdx==2);
    end
    
    % dispersion
    omega = [omega(:,kmIdx==1); omega(:,kmIdx==2); omega(:,kmIdx==3)];
    % exchange matrices
    Sab   = cat(3,mmat(Sab(:,:,:,kmIdx==1),K1), mmat(Sab(:,:,:,kmIdx==2),K2), ...
        mmat(Sab(:,:,:,kmIdx==3),conj(K1)));
    
    hkl   = hkl(:,kmIdx==2);
    nHkl0 = nHkl0/3;
else
    helical = false;
end

if ~param.notwin
    if nTwin > 1
        omega = mat2cell(omega,size(omega,1),repmat(nHkl0,[1 nTwin]));
    end
    % Rotate the calculated correlation function into the twin coordinate
    % system using rotC
    SabAll = cell(1,nTwin);
    for ii = 1:nTwin
        % select the ii-th twin from the Q points
        idx    = (1:nHkl0) + (ii-1)*nHkl0;
        % select correlation function of twin ii
        SabT   = Sab(:,:,:,idx);
        if param.sortMode
            % sort the spin wave modes
            [omega{ii}, SabT] = sortmode(omega{ii},reshape(SabT,9,size(SabT,3),[]));
            SabT = reshape(SabT,3,3,size(SabT,2),[]);
        end
        % size of the correlation function matrix
        sSabT  = size(SabT);
        % convert the matrix into cell of 3x3 matrices
        SabT   = reshape(SabT,3,3,[]);
        % select the rotation matrix of twin ii
        rotC   = obj.twin.rotc(:,:,ii);
        % rotate correlation function using arrayfun
        SabRot = arrayfun(@(idx)(rotC*SabT(:,:,idx)*(rotC')),1:size(SabT,3),'UniformOutput',false);
        SabRot = cat(3,SabRot{:});
        % resize back the correlation matrix
        SabAll{ii} = reshape(SabRot,sSabT);
    end
    Sab = SabAll;
    
    if nTwin == 1
        Sab = Sab{1};
    end
else
    if param.sortMode
        % sort the spin wave modes
        [omega, Sab] = sortmode(omega,reshape(Sab,9,size(Sab,3),[]));
        Sab          = reshape(Sab,3,3,size(Sab,2),[]);
    end
end

% Creates output structure with the calculated values.
spectra.omega    = omega;
spectra.Sab      = Sab;
spectra.hkl      = obj.unit.qmat\hkl(:,1:nHkl0);
spectra.hklA     = hklA;
spectra.incomm   = incomm;
spectra.helical  = helical;
spectra.norm     = false;
spectra.nformula = double(obj.unit.nformula);

% Save different intermediate results.
if param.saveV
    if param.notwin
        spectra.V = Vsave;
    else
        spectra.V = mat2cell(Vsave, size(Vsave, 1), size(Vsave, 1), repmat(nHkl0, [1 nTwin]));
    end
end
if param.saveH
    spectra.H = Hsave;
end
if param.saveSabp && incomm
    spectra.Sabp = Sabp;
    spectra.omegap = omegap;
end

% save the important parameters
spectra.param.notwin    = param.notwin;
spectra.param.sortMode  = param.sortMode;
spectra.param.tol       = param.tol;
spectra.param.omega_tol = param.omega_tol;
spectra.param.hermit    = param.hermit;
spectra.title           = param.title;
spectra.gtensor         = param.gtensor;

if ~param.fitmode
    spectra.dateend = datestr(now);
    spectra.obj = copy(obj);
end

if ~param.gtensor && any(obj.single_ion.g)
    warning('spinw:spinwave:NonZerogTensor',['The SpinW model defines a '...
        'g-tensor that is not included in the calculation. Anisotropic '...
        'g-tensor values cannot be applied afterwards as they change relative'...
        'spin wave intensities!'])
end

% issue eigorth warning
if orthWarn0
    warning('spinw:spinwave:NoOrth','Eigenvectors of defective eigenvalues cannot be orthogonalised at some q-point!');
end

if strcmp(singWarn,'MATLAB:nearlySingularMatrix')
    lineLink = 'line 846';
    if feature('HotLinks')
        lineLink = ['<a href="matlab:opentoline([''' sw_rootdir 'swfiles' filesep '@spinw' filesep 'spinwave.m''' '],846,0)">' lineLink '</a>'];
    end
    warning('spinw:spinwave:nearlySingularMatrix',['Matrix is close '...
        'to singular or badly scaled. Results may be inaccurate.\n> In spinw/spinwave (' lineLink ')']);
    %fprintf(repmat('\b',[1 30]));
end

end
