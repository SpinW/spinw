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
% `'cmplxBase'`
% : If `true`, we use a local coordinate system fixed by the
%   complex magnetisation vectors:
%   $\begin{align}  e_1 &= \Im(\hat{M})\\
%                   e_3 &= Re(\hat{M})\\
%                   e_2 &= e_3\times e_1
%    \end{align}$
%   If `false`, we use a coordinate system fixed to the moments:
%   $\begin{align}  e_3 \parallel S_i\\
%                   e_2 &= \S_i \times [1, 0, 0]\\
%                   e_1 &= e_2 \times e_3
%   \end{align}$
%   Except if $S_i \parallel [1, 0, 0], e_2 = [0, 0, 1]$. The default is
%  `false`.
%
% `'gtensor'`
% : If true, the g-tensor will be included in the spin-spin correlation
%   function. Including anisotropic g-tensor or different
%   g-tensor for different ions is only possible here. Including a simple
%   isotropic g-tensor is possible afterwards using the [sw_instrument]
%   function.
%
% `'neutron_output'`
% : If `true`, will output only `Sperp`, the S(q,w) component perpendicular
%   to Q that is measured by neutron scattering, and will *not* output the
%   full Sab tensor. (Usually sw_neutron is used to calculate `Sperp`.)
%   Default value is `false`.
%
% `'fitmode'`
% : If `true`, function is optimized for multiple consecutive calls (e.g.
%   the output spectrum won't contain the copy of `obj`), default is
%   `false`.
%
% `'fastmode'`
% : If `true`, will set `'neutron_output', true`, `'fitmode', true`,
%   `'sortMode', false`, and will only output intensity for positive energy
%   (neutron energy loss) modes. Default value is `false`.
%
% `'notwin'`
% : If `true`, the spectra of the twins won't be calculated. Default is
%   `false`.
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
%   $S'(k,\omega)$ is saved. For incommensurate structures only. Default
%   value is `false`.
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
%   * `Sperp`   The component of `Sab` perpendicular to $Q$, which neutron
%               scattering measures. This is outputed *instead* of `Sab`
%               if the `'neutron_output', true` is specified.
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
%   * `hkl`     Contains the input $Q$ values, dimensions are $[3\times n_{Q}]$.
%   * `hklA`    Same $Q$ values, but in $\\ang^{-1}$ unit, in the
%               lab coordinate system, dimensins are $[3\times n_{Q}]$.
%   * `formfact`Logical value, whether the form factor has been included in
%               the spin-spin correlation function.
%   * `incomm`  Logical value, tells whether the calculated spectra is
%               incommensurate or not.
%   * `helical` Logical value, whether the magnetic structure is a helix
%               i.e. whether 2*k is non-integer.
%   * `norm`    Logical value, is always false.
%   * `nformula`Number of formula units in the unit cell that have been
%               used to scale Sab, as given in spinw.unit.nformula.
%   * `param`   Struct containing input parameters, each corresponds to the
%               input parameter of the same name:
%               * `notwin`
%               * `sortMode`
%               * `tol`
%               * `omega_tol`
%               * `hermit`
%   * `title`   Character array, the title for the output spinwave, default
%               is 'Numerical LSWT spectrum'
%   * `gtensor` Logical value, whether a g-tensor has been included in the
%               calculation.
%   * `obj`     The copy (clone) of the input `obj`, see [spinw.copy].
%   * `datestart`Character array, start date and time of the calculation
%   * `dateend` Character array, end date and time of the calculation
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

if isa(hkl, 'sym') && ~obj.symbolic
    obj.symbolic(true);
    setnosym = onCleanup(@()obj.symbolic(false));
end

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
        spectra = obj.spinwavesym(varargin{:},'hkl',hkl(:));
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

inpForm.fname  = [inpForm.fname  {'neutron_output' 'fastmode'}];
inpForm.defval = [inpForm.defval {false             false }];
inpForm.size   = [inpForm.size   {[1 1]             [1 1] }];

param = sw_readparam(inpForm, varargin{:});

if param.fastmode
    param.neutron_output = true;
    param.fitmode = true;
    param.sortMode = false;
    if any(useMex) && (param.saveV || param.saveH || param.saveSabp)
        warning('spinw:spinwave:fastmodewithsave', ...
                ['You have set both "usemex" and "fastmode" and also ' ...
                 'requested that S, V or H is saved, but mex files do not ' ...
                 'support saving intermediate matrices. So in this case ' ...
                 'the mex files will *not* be used. Set all "save*" ' ...
                 'options to "false" to use mex files.']);
    end
end

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
if incomm && prod(nExt) > 1
    warning('spinw:spinwave:IncommKinSupercell', ...
            ['The results for an incommensurate modulation in a ' ...
             'supercell have not been scientifically validated.']);
end
if ~incomm && param.saveSabp
    warning('spinw:spinwave:CommensurateSabp', ['The dynamical structure '...
            'factor in the rotating frame has been requested, but the ', ...
            'structure is commensurate so this will have no effect.']);
end

% If only one hkl value, convert to column vector
if numel(hkl) == 3
    hkl = hkl(:);
end

% Transform the momentum values to the new lattice coordinate system
hkl = obj.unit.qmat*hkl;

% Calculates momentum transfer in A^-1 units.
hklA = 2*pi*(hkl'/obj.basisvector)';

% Check for 2*km (if 2*km is integer, the magnetic structure is not a true helix)
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

    if param.saveSabp
        Sabp = [];
        omegap = [];
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
        if param.neutron_output
            hklAf{tt} = [hklA hklA hklA];
        end
    end
    nHkl  = nHkl*3;
    nHkl0 = nHkl0*3;
else
    hkl0   = hkl;
    hklExt = hkl;
    helical = false;
end

hkl    = cell2mat(hkl);
hkl0   = cell2mat(hkl0);
hklExt = cell2mat(hklExt);

if param.neutron_output
    if incomm
        hklAf = cell2mat(hklAf);
    elseif ~param.notwin
        hklAf = repmat(hklA, [1 nTwin]);
    else
        hklAf = hklA;
    end
    % Normalized scattering wavevector in xyz coordinate system.
    hklAf = bsxfun(@rdivide, hklAf, sqrt(sum(hklAf.^2, 1)));
else
    hklAf = [];
end
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
    % e3 || Si
    e3 = bsxfun(@rdivide,M0,S0);
    % e2 = Si x [1,0,0], if Si || [1,0,0] --> e2 = [0,0,1]
    e2  = [zeros(1,nMagExt); e3(3,:); -e3(2,:)];
    e2(3,~any(abs(e2)>1e-10)) = 1;
    e2  = bsxfun(@rdivide,e2,sqrt(sum(e2.^2,1)));
    % e1 = e2 x e3
    e1  = cross(e2,e3);
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
else
    bqdR = [];
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

z1 = zed;
if param.gtensor

    gtensor = SI.g;

    if incomm
        % keep the rotation invariant part of g-tensor
        nx  = [0 -n(3) n(2);n(3) 0 -n(1);-n(2) n(1) 0];
        nxn = n'*n;
        m1  = eye(3);
        gtensor = 1/2*gtensor - 1/2*mmat(mmat(nx,gtensor),nx) + 1/2*mmat(mmat(nxn-m1,gtensor),nxn) + 1/2*mmat(mmat(nxn,gtensor),2*nxn-m1);
    end
    for i1 = 1:size(gtensor, 3)
        z1(:,i1) = gtensor(:,:,i1) * zed(:,i1);
    end
end

hklIdx = [floor(((1:nSlice)-1)/nSlice*nHkl)+1 nHkl+1];

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
    FF = [];
end

if (~isstring(useMex) && ~ischar(useMex) && useMex) || strcmp(useMex, 'auto')
    % For large unit cells, the original mex code is much better because it parallelizes the individual
    % eigen/cholesky decompositions operations (chol / eig) in addition to parallelising over Q-points
    % The new code uses Eigen for these operations which is strictly serial so parallises better over Q-points
    % but will be super slow for large nMagExt
    if nMagExt > pref.nspinlarge % This threshold needs to be explored more
        useMex = 'old';
    else
        useMex = 'new';
    end
end

use_swloop = any(useMex) && strcmp(useMex, 'new') && ~param.saveH && ~param.saveV && ~param.saveSabp;

if ~use_swloop
    if param.fastmode
        nMode = nMagExt;
    else
        nMode = 2 * nMagExt;
    end

    % Empty omega dispersion of all spin wave modes, size: 2*nMagExt x nHkl.
    omega = zeros(2*nMagExt, nHkl);
    if param.neutron_output
        Sperp = zeros(nMode, nHkl);
    else
        SabFull = zeros(3,3,nMode,nHkl);
    end

    if incomm
        % resize matrices due to the incommensurability (k-km,k,k+km) multiplicity
        kmIdx = repmat(sort(repmat([1 2 3],1,nHkl0/3)),1,nTwin);
        % Rodrigues' rotation formula.
        nx  = [0 -n(3) n(2); n(3) 0 -n(1); -n(2) n(1) 0];
        nxn = n'*n;
        K1 = 1/2*(eye(3) - nxn - 1i*nx);
        K2 = nxn;
        m1  = eye(3);
    end
end

if use_swloop
    pars = struct('hermit', param.hermit, 'omega_tol', param.omega_tol, 'formfact', param.formfact, ...
        'incomm', incomm, 'helical', helical, 'nTwin', nTwin, 'bq', any(bq), 'field', any(SI.field), ...
        'nThreads', pref.nthread, 'n', n, 'rotc', obj.twin.rotc, 'fastmode', param.fastmode, ...
        'neutron_output', param.neutron_output, 'hklA', hklAf, 'nformula', obj.unit.nformula, ...
        'nCell', prod(nExt));
    ham_diag = diag(accumarray([idxA2; idxD2], 2*[A20 D20], [1 1]*2*nMagExt));
    idxAll = [idxA1; idxB; idxD1]; ABCD = [AD0 2*BC0 conj(AD0)];
    bqABCD = []; bq_ham_d = []; idxBq = []; ham_MF = {};
    if pars.bq
        bqABCD = [bqA0 conj(bqA0) 2*bqB0];
        idxBq = [idxbqA; idxbqA2; idxbqB];
        bq_ham_d = diag(accumarray([idxbqC; idxbqC2; idxbqD], [bqC bqC 2*bqD], [1 1]*2*nMagExt));
        assert(isreal(bqABCD), 'Internal logical error');
    end
    if pars.field
        for ii = 1:nTwin
            ham_MF{ii} = accumarray(idxMF, MF(:,:,ii), [1 1]*2*nMagExt);
        end
    end
    try
        [omega, Sab, warn1, orthWarn0] = swloop(pars, hklExt, ...
            ABCD, idxAll, ham_diag, dR, RR, S0, z1, FF, bqdR, bqABCD, idxBq, bq_ham_d, ham_MF);
    catch err
        if ~isempty(strfind(err.message, 'notposdef'))
            error('chol_omp:notposdef', 'Hamiltonian is not positive definite');
        elseif ~isempty(strfind(err.message, 'Eigensolver'))
            error('swloop:notconverge', 'Could not determine eigenvalues of Hamiltonian');
        else
            rethrow(err);
        end
    end
    if param.hermit && sum(abs(imag(omega(:)))) < 1e-5
        omega = real(omega);
    end
    if param.neutron_output
        Sperp = Sab;
        clear Sab;
    end
    if incomm
        nHkl0 = nHkl0/3;
        nHklT = nHkl / nTwin;
        kmIdx = cell2mat(arrayfun(@(x) (x+nHkl0):(x+2*nHkl0-1), [0:(nTwin-1)]*nHklT + 1, 'UniformOutput', false));
        hkl = hkl(:, kmIdx);
    end
else
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

        if any(useMex) && nHklMEM>1
            % use mex files to speed up the calculation
            % mex file will return an error if the matrix is not positive definite.
            [Ksq, invK] = chol_omp(ham,'Colpa','tol',param.omega_tol);
            [V, omega(:,hklIdxMEM)] = eig_omp(Ksq,'sort','descend');
            % the inverse of the para-unitary transformation V
            if param.fastmode
                % Only transform the positive energy modes (first half of V)
                for ii = 1:nMagExt
                    V(:,ii,:) = bsxfun(@times, squeeze(V(:,ii,:)), sqrt(omega(ii,hklIdxMEM)));
                end
                V = sw_mtimesx(invK, V(:,1:nMagExt,:));
            else
                for ii = 1:nHklMEM
                    V(:,:,ii) = V(:,:,ii)*diag(sqrt(gCommd.*omega(:,hklIdxMEM(ii))));
                end
                % V = bsxfun(@times, invK, V);
                V = sw_mtimesx(invK, V);
            end
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

                Ksq = K*gComm*K';
                Ksq = 1/2*(Ksq+Ksq');
                % Hermitian Ksq will give orthogonal eigenvectors
                [U, D] = eig(Ksq);
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
        if any(useMex)
            gham = sw_mtimesx(gComm, ham);
        else
            gham = mmat(gComm,ham);
        end

        [V, D, orthWarn] = eigorth(gham,param.omega_tol,useMex);

        orthWarn0 = orthWarn || orthWarn0;

        if param.fastmode
            % Only transform the positive energy modes (first half of V)
            for ii = 1:nMagExt
                V(:,ii,:) = bsxfun(@times, V(:,ii,:), sqrt(1 ./ sum(bsxfun(@times,gCommd,conj(V(:,ii,:)).*V(:,ii,:)))));
            end
        else
            for ii = 1:nHklMEM
                % multiplication with g removed to get negative and positive
                % energies as well
                omega(:,hklIdxMEM(ii)) = D(:,ii);
                M              = diag(gComm*V(:,:,ii)'*gComm*V(:,:,ii));
                V(:,:,ii)      = V(:,:,ii)*diag(sqrt(1./M));
            end
        end
    end

    if param.saveV
        Vsave(:,:,hklIdxMEM) = V;
    end
    if param.saveH
        Hsave(:,:,hklIdxMEM) = ham;
    end

    if param.fastmode
        V = V(:,1:nMagExt,:);
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

    ExpFL =      repmat(permute(ExpF,[1 4 5 2 3]),[3 3 nMode 2]);
    % conj transpose of ExpFL
    ExpFR = conj(permute(ExpFL,[1 2 4 3 5]));

    zeda = repmat(permute([zed conj(zed)],[1 3 4 2]),[1 3 nMode 1 nHklMEM]);
    % conj transpose of zeda
    zedb = conj(permute(zeda,[2 1 4 3 5]));

    % calculate magnetic structure factor using the hklExt0 Q-values
    % since the S(Q+/-k,omega) correlation functions also belong to the
    % F(Q)^2 form factor

    if param.formfact
        % include the form factor in the z^alpha, z^beta matrices
        zeda = zeda.*repmat(permute(FF(:,hklIdxMEM),[3 4 5 1 2]),[3 3 nMode 2 1]);
        zedb = zedb.*repmat(permute(FF(:,hklIdxMEM),[3 4 1 5 2]),[3 3 2 nMode 1]);
    end

    if param.gtensor
        % include the g-tensor
        zeda = mmat(repmat(permute(gtensor,[1 2 4 3]),[1 1 1 2]),zeda);
        zedb = mmat(zedb,repmat(gtensor,[1 1 2]));
    end


    % Dynamical structure factor from S^alpha^beta(k) correlation function.
    % Sab(alpha,beta,iMode,iHkl), size: 3 x 3 x 2*nMagExt x nHkl.
    % Normalizes the intensity to single unit cell.
    Sab = reshape(sum(zeda.*ExpFL.*VExtL,4),[3 3 nMode nHklMEM]) .* ...
          reshape(sum(zedb.*ExpFR.*VExtR,3),[3 3 nMode nHklMEM]) / prod(nExt);

    if incomm
        if helical
            % integrating out the arbitrary initial phase of the helix
            Sab = 1/2*Sab - 1/2*mmat(mmat(nx,Sab),nx) + 1/2*mmat(mmat(nxn-m1,Sab),nxn) + 1/2*mmat(mmat(nxn,Sab),2*nxn-m1);
        end
        kmIdxMEM = kmIdx(hklIdxMEM);

        % Save the structure factor in the rotating frame
        if param.saveSabp
            Sabp = cat(4, Sabp, Sab(:,:,:,kmIdxMEM==2));
            omegap = cat(2, omega(:,kmIdxMEM==2));
        end

        % Rotate Sab back into lab frame
        Sab(:,:,:,kmIdxMEM==1) = mmat(Sab(:,:,:,kmIdxMEM==1), K1);
        Sab(:,:,:,kmIdxMEM==2) = mmat(Sab(:,:,:,kmIdxMEM==2), K2);
        Sab(:,:,:,kmIdxMEM==3) = mmat(Sab(:,:,:,kmIdxMEM==3), conj(K1));
    end

    if ~param.notwin
        % Rotate the calculated correlation function into the twin coordinate system using rotC
        for ii = 1:nTwin
            % select the ii-th twin from the Q points
            idx = find((hklIdxMEM <= (nHkl0 * ii)) .* (hklIdxMEM > (nHkl0 * (ii-1))));
            if ~isempty(idx)
                % convert the matrix into cell of 3x3 matrices
                SabT   = reshape(Sab(:,:,:,idx),3,3,[]);
                % select the rotation matrix of twin ii
                rotC   = obj.twin.rotc(:,:,ii);
                % rotate correlation function using arrayfun
                SabRot = arrayfun(@(idx)(rotC*SabT(:,:,idx)*(rotC')),1:size(SabT,3),'UniformOutput',false);
                SabRot = cat(3,SabRot{:});
                % resize back the correlation matrix
                Sab(:,:,:,idx) = reshape(SabRot, [3 3 nMode numel(idx)]);
            end
        end
    end

    if param.neutron_output
        if obj.unit.nformula > 0
            Sab = Sab/double(obj.unit.nformula);
        end

        % get symmetric component of Sab only
        Sab = (Sab + permute(Sab,[2 1 3 4]))/2;

        % Normalized scattering wavevector in xyz coordinate system.
        hklAN = hklAf(:, hklIdxMEM);

        % avoid NaN for Q=0
        NaNidx = find(any(isnan(hklAN)));
        for kk = 1:numel(NaNidx)
            if NaNidx(kk) < size(hklAN,2)
                hklAN(:,NaNidx(kk)) = hklAN(:,NaNidx(kk)+1);
            else
                hklAN(:,NaNidx(kk)) = [1;0;0];
            end
        end

        hkla = repmat(permute(hklAN,[1 3 2]),[1 3 1]);
        hklb = repmat(permute(hklAN,[3 1 2]),[3 1 1]);

        % Perpendicular part of the scattering wavevector.
        qPerp = repmat(eye(3),[1 1 numel(hklIdxMEM)])- hkla.*hklb;
        qPerp = repmat(permute(qPerp,[1 2 4 3]),[1 1 nMode 1]);

        % Dynamical structure factor for neutron scattering
        % Sperp: nMode x nHkl.
        Sperp(:,hklIdxMEM) = permute(sumn(qPerp.*Sab,[1 2]),[3 4 1 2]);
    else
        SabFull(:,:,:,hklIdxMEM) = Sab;
    end

    sw_timeit(jj/nSlice*100,0,param.tid);
  end
  if param.fastmode
    omega = omega(1:nMagExt,:);
  end
  if ~param.neutron_output
    Sab = SabFull;
  end
end

[~,singWarn] = lastwarn;
% restore warning for singular matrix
warning(singWarn0.state,'MATLAB:nearlySingularMatrix');

% If number of formula units are given per cell normalize to formula
% unit
if obj.unit.nformula > 0 && ~param.neutron_output
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

if incomm && ~use_swloop
    % Rearrange the Sab and omega matrices
    omega = [omega(:,kmIdx==1); omega(:,kmIdx==2); omega(:,kmIdx==3)];
    if ~param.neutron_output
        Sab = cat(3, Sab(:,:,:,kmIdx==1), Sab(:,:,:,kmIdx==2), Sab(:,:,:,kmIdx==3));
    else
        Sperp = [Sperp(:,kmIdx==1); Sperp(:,kmIdx==2); Sperp(:,kmIdx==3)];
    end
    hkl   = hkl(:,kmIdx==2);
    nHkl0 = nHkl0/3;
end

if ~param.notwin && nTwin > 1
    omega = mat2cell(omega, size(omega, 1), repmat(nHkl0, [1 nTwin]));
    if param.neutron_output
        Sperp = mat2cell(Sperp, size(Sperp, 1), repmat(nHkl0, [1 nTwin]));
    else
        Sab = squeeze(mat2cell(Sab, 3, 3, size(Sab, 3), repmat(nHkl0, [1 nTwin])))';
    end
end

if param.sortMode && ~param.neutron_output
    if ~param.notwin
        for ii = 1:nTwin
            % sort the spin wave modes
            [omega{ii}, Sab{ii}] = sortmode(omega{ii},reshape(Sab{ii},9,size(Sab{ii},3),[]));
            Sab{ii} = reshape(Sab{ii},3,3,size(Sab{ii},2),[]);
        end
    else
        % sort the spin wave modes
        [omega, Sab] = sortmode(omega,reshape(Sab,9,size(Sab,3),[]));
        Sab = reshape(Sab,3,3,size(Sab,2),[]);
    end
end

% Creates output structure with the calculated values.
spectra.omega    = omega;
if param.neutron_output
    spectra.Sperp = Sperp;
else
    spectra.Sab  = Sab;
end
spectra.hkl      = obj.unit.qmat\hkl(:,1:nHkl0);
spectra.hklA     = hklA;
spectra.incomm   = incomm;
spectra.helical  = helical;
spectra.norm     = false;
spectra.nformula = int32(obj.unit.nformula);

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

if param.fitmode
    % Copies only fields needed by downstream functions (sw_egrid, sw_neutron, sw_plotspec)
    spectra.obj = struct('single_ion', obj.single_ion, 'twin', obj.twin, ...
        'unit', obj.unit, 'basisvector', obj.basisvector, 'nmagext', obj.nmagext);
else
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
    warning('spinw:spinwave:nearlySingularMatrix',['Matrix is close '...
            'to singular or badly scaled. Results may be inaccurate.']);
end

end
