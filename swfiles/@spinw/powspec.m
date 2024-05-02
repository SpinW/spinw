function spectra = powspec(obj, hklA, varargin)
% calculates powder averaged spin wave spectra
%
% ### Syntax
%
% `spectra = powspec(obj,QA)`
%
% `spectra = powspec(___,Name,Value)`
%
% ### Description
%
% `spectra = powspec(obj,QA)` calculates powder averaged spin wave spectrum
% by averaging over spheres with different radius around origin in
% reciprocal space. This way the spin wave spectrum of polycrystalline
% samples can be calculated. This method is not efficient for low
% dimensional (2D, 1D) magnetic lattices. To speed up the calculation with
% mex files use the `swpref.setpref('usemex',true)` option. 
%
% `spectra = powspec(___,Value,Name)` specifies additional parameters for
% the calculation. For example the function can calculate powder average of
% arbitrary spectral function, if it is specified using the `specfun`
% option. 
%
% ### Example
%
% Using only a few lines of code one can calculate the powder spectrum of
% the triangular lattice antiferromagnet ($S=1$, $J=1$) between $Q=0$ and 3
% \\ang$^{-1}$ (the lattice parameter is 3 \\ang).
%
% ```
% >>tri = sw_model('triAF',1);
% >>E = linspace(0,4,100);
% >>Q = linspace(0,4,300);
% >>triSpec = tri.powspec(Q,'Evect',E,'nRand',1e3);
% >>sw_plotspec(triSpec);
% >>snapnow
% ```
%
% ### Input arguments
%
% `obj`
% : [spinw] object.
%
% `QA`
% : Vector containing the $Q$ values in units of the inverse of the length
% unit (see [spinw.unit]) with default unit being \\ang$^{-1}$. The
% value are stored in a row vector with $n_Q$ elements.
%
% ### Name-Value Pair Arguments
%
% `specfun`
% : Function handle of a solver. Default value is `@spinwave`. It is
%   currently tested with two functions:
%
%   * `spinw.spinwave` 	Powder average spin wave spectrum.
%   * `spinw.scga`      Powder averaged diffuse scattering spectrum.
%
% `nRand`
% : Number of random orientations per `QA` value, default value is 100.
%
% `Evect`
% : Row vector, defines the center/edge of the energy bins of the
%   calculated output, number of elements is $n_E$. The energy units are
%   defined by the `spinw.unit.kB` property. Default value is an edge bin
%   `linspace(0,1.1,101)`.
%
% `binType`
% : String, determines the type of bin, possible options:
%   * `'cbin'`    Center bin, the center of each energy bin is given.
%   * `'ebin'`    Edge bin, the edges of each bin is given.
%
%   Default value is `'ebin'`.
%
% `'T'`
% : Temperature to calculate the Bose factor in units
%   depending on the Boltzmann constant. Default value taken from
%   `obj.single_ion.T` value.
%
% `'title'`
% : Gives a title to the output of the simulation.
%
% `'extrap'`
% : If true, arbitrary additional parameters are passed over to
%   the spectrum calculation function.
%
% `'fibo'`
% : If true, instead of random sampling of the unit sphere the Fibonacci
%   numerical integration is implemented as described in
%   [J. Phys. A: Math. Gen. 37 (2004)
%   11591](http://iopscience.iop.org/article/10.1088/0305-4470/37/48/005/meta).
%   The number of points on the sphere is given by the largest
%   Fibonacci number below `nRand`. Default value is false.
%
% `'imagChk'`
% : Checks that the imaginary part of the spin wave dispersion is
%   smaller than the energy bin size. Default value is true.
%
% `'component'`
% : See [sw_egrid] for the description of this parameter.
%
% The function also accepts all parameters of [spinw.spinwave] with the
% most important parameters are:
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
% `dE`
% : Energy resolution (FWHM) can be function, or a numeric matrix that
%   has length 1 or the number of energy bin centers.
%
% `'neutron_output'`
% : If `true`, the spinwave will output only `Sperp`, the S(q,w) component
%   perpendicular to Q that is measured by neutron scattering, and will
%   *not* output the  full Sab tensor. (Usually sw_neutron is used to 
%   calculate `Sperp`.) Default value is `false`.
%
% `'fastmode'`
% : If `true`, will set `'neutron_output', true`, `'fitmode', true`,
%   `'sortMode', false`, and will only output intensity for positive energy
%   (neutron energy loss) modes. Default value is `false`.
%
% The function accepts some parameters of [spinw.scga] with the most important
% parameters are:
%
% `'nInt'`
% : Number of $Q$ points where the Brillouin zone is sampled for the
%   integration.
%
% ### Output Arguments
%
% `spectra`
% : structure with the following fields:
%
%   * `swConv` The spectra convoluted with the dispersion. The center
%     of the energy bins are stored in `spectra.Evect`. Dimensions are
%     $[n_E\times n_Q]$.
%   * `hklA` Same $Q$ values as the input `hklA`.
%   * `Evect` Contains the bins (edge values of the bins) of the energy transfer
%     values, dimensions are $[1\times n_E+1]$.
%   * `param` Contains all the input parameters.
%   * `obj` The clone of the input `obj` object, see [spinw.copy].
%
% ### See also
%
% [spinw] \| [spinw.spinwave] \| [spinw.optmagstr] \| [sw_egrid]
%

% help when executed without argument
if nargin==1
    swhelp spinw.powspec
    return
end

pref = swpref;

hklA = hklA(:)';
T0 = obj.single_ion.T;

title0 = 'Powder LSWT spectrum';
tid0   = pref.tid;

inpForm.fname  = {'nRand' 'Evect'    'T'   'formfact' 'formfactfun' 'tid' 'nInt'};
inpForm.defval = {100     zeros(1,0) T0    false      @sw_mff       tid0  1e3   };
inpForm.size   = {[1 1]   [1 -1]     [1 1] [1 -2]     [1 1]         [1 1] [1 1] };
inpForm.soft   = {false   false      false false      false         false false   };

inpForm.fname  = [inpForm.fname  {'hermit' 'gtensor' 'title' 'specfun' 'imagChk'}];
inpForm.defval = [inpForm.defval {true     false     title0  @spinwave  true    }];
inpForm.size   = [inpForm.size   {[1 1]    [1 1]     [1 -3]  [1 1]      [1 1]   }];
inpForm.soft   = [inpForm.soft   {false    false     false   false      false}];

inpForm.fname  = [inpForm.fname  {'extrap' 'fibo' 'optmem' 'binType' 'component'}];
inpForm.defval = [inpForm.defval {false    false  0        'ebin'    'Sperp'    }];
inpForm.size   = [inpForm.size   {[1 1]    [1 1]  [1 1]    [1 -4]     [1 -5]    }];
inpForm.soft   = [inpForm.soft   {false    false  false    false      false}];

inpForm.fname  = [inpForm.fname  {'fid' , 'dE'}];
inpForm.defval = [inpForm.defval {-1,      []}];
inpForm.size   = [inpForm.size   {[1 1],   [-6, -7]}];
inpForm.soft   = [inpForm.soft   {false    true}];

inpForm.fname  = [inpForm.fname  {'neutron_output' 'fastmode'}];
inpForm.defval = [inpForm.defval {false             false }];
inpForm.size   = [inpForm.size   {[1 1]             [1 1] }];
inpForm.soft   = [inpForm.soft   {false             false}];


param  = sw_readparam(inpForm, varargin{:});

if param.fid == -1
    fid = pref.fid;
else
    fid = param.fid;
end

% list of supported functions:
%   0:  unknown
%   1:  @spinwave
%   2:  @scga
funList = {@spinwave @scga};
funIdx  = [find(cellfun(@(C)isequal(C,param.specfun),funList)) 0];
funIdx  = funIdx(1);

if isempty(param.Evect) && funIdx == 1
    error('spinw:powspec:WrongOption','Energy bin vector is missing, use ''Evect'' option!');
end

% number of bins along energy
switch param.binType
    case 'cbin'
        nE      = numel(param.Evect);
    case 'ebin'
        nE      = numel(param.Evect) - 1;
end

nQ      = numel(hklA);

fprintf0(fid,'Calculating powder spectra...\n');

% message for magnetic form factor calculation
yesNo = {'No' 'The'};
fprintf0(fid,[yesNo{param.formfact+1} ' magnetic form factor is'...
    ' included in the calculated structure factor.\n']);
% message for g-tensor calculation
fprintf0(fid,[yesNo{param.gtensor+1} ' g-tensor is included in the '...
    'calculated structure factor.\n']);

sw_timeit(0,1,param.tid,'Powder spectrum calculation');

nk = param.nRand;
if param.fibo
    % apply the Fibonacci numerical integration on a sphere
    % according to J. Phys. A: Math. Gen. 37 (2004) 11591
    % create QF points on the unit sphere
    
    [F,F1] = fibonacci(param.nRand);
    param.nRand = F;
    
    QF = zeros(3,F);
    
    j = 0:(F-1);
    QF(3,:) = j*2/F-1;
    
    theta = asin(QF(3,:));
    phi   = 2*pi*F1/F*j;
    
    QF(1,:) = cos(theta).*sin(phi);
    QF(2,:) = cos(theta).*cos(phi);
    
    nk = F;
end

% lambda value for SCGA, empty will make integration in first loop
specQ.lambda = [];

hkl = zeros(3, nk * nQ);
for ii = 1:nQ
    if param.fibo
        Q = QF*hklA(ii);
    else
        rQ  = randn(3,param.nRand);
        Q   = bsxfun(@rdivide,rQ,sqrt(sum(rQ.^2)))*hklA(ii);
    end
    i0 = (ii-1) * nk + 1;
    hkl(:, i0:(i0+nk-1)) = (Q'*obj.basisvector)'/2/pi;
end

switch funIdx
    case 0
        % general function call allow arbitrary additional parameters to
        % pass to the spectral calculation function
        warnState = warning('off','sw_readparam:UnreadInput');
        specQ = param.specfun(obj,hkl,varargin{:});
        warning(warnState);
        specQ = sw_neutron(specQ,'pol',false);
    case 1
        % @spinwave
        specQ = spinwave(obj,hkl,struct('fitmode',true,'notwin',true,...
            'Hermit',param.hermit,'formfact',param.formfact,...
            'formfactfun',param.formfactfun,'gtensor',param.gtensor,...
            'optmem',param.optmem,'tid',param.tid,'fid',0, ...
            'neutron_output', param.neutron_output, 'fastmode', ...
            param.fastmode),'noCheck');
    case 2
        % @scga
        specQ = scga(obj,hkl,struct('fitmode',true,'formfact',param.formfact,...
            'formfactfun',param.formfactfun,'gtensor',param.gtensor,...
            'fid',0,'lambda',specQ.lambda,'nInt',param.nInt,'T',param.T,...
            'plot',false),'noCheck');
        specQ = sw_neutron(specQ,'pol',false);
end
specQ.obj = obj;
% use edge grid by default
specQ = sw_egrid(specQ,struct('Evect',param.Evect,'T',param.T,'binType',param.binType,...
    'imagChk',param.imagChk,'component',param.component, 'dE', param.dE),'noCheck');
powSpec = squeeze(sum(reshape(specQ.swConv, max(1, nE), nk, nQ), 2) / param.nRand);

fprintf0(fid,'Calculation finished.\n');

% save different field into spectra
spectra.swConv    = powSpec;
spectra.hklA      = hklA;
spectra.component = param.component;
spectra.nRand    = param.nRand;
spectra.T        = param.T;
spectra.obj      = copy(obj);
spectra.norm     = false;
spectra.formfact = specQ.formfact;
spectra.gtensor  = specQ.gtensor;
spectra.date     = datestr(now);
spectra.title    = param.title;
% save all input parameters of spinwave into spectra
spectra.param    = specQ.param;

% some spectral function dependent parameters
switch funIdx
    case 0
        spectra.Evect    = specQ.Evect;
    case 1
        spectra.Evect    = specQ.Evect;
        spectra.incomm   = specQ.incomm;
        spectra.helical  = specQ.helical;
    case 2
        spectra.lambda   = specQ.lambda;
end

end

function [F,F1] = fibonacci(Fmax)
% returns the last two Fibonacci number smaller or equal to the
% given number
%
% [Flast Fprev] = fibonacci(Fmax)
%

num = [0 0 1];

while num(end)<Fmax
    num(end+1) = sum(num(end+[-1 0])); %#ok<AGROW>
end

if num(end) == Fmax
    F = num(end);
    F1 = num(end-1);
else
    F = num(end-1);
    F1 = num(end-2);
end
end
