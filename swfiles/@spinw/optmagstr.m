function optm = optmagstr(obj, varargin)
% general magnetic structure optimizer
% 
% ### Syntax
% 
% `optm = optmagstr(obj,Name,Value)`
% 
% ### Description
% 
% `optm = optmagstr(obj,Name,Value)` is a general magnetic structure
% optimizer that as the name suggests is general however usually less
% efficient than [spinw.optmagk] or [spinw.optmagsteep]. However this
% function enables the usage of constraint functions to improve the
% optimization. This function is most useful if there are 1-2 parameters
% that have to be optimized, such as a canting angle of the spins in
% a magnetic field. To optimize large numbers of spin angles
% [spinw.optmagsteep] might be faster. Only obj.mag_str.nExt is used from
% an already initialised magnetic structure, initial k and S are determined
% from the optimisation function parameters. If a magnetic structure has
% not been initialised in obj, nExt = [1 1 1] is used.
% 
% ### Examples
% 
% The example determines the propagation vector of the ground state of the
% triangular lattice antiferromagnet. The magnetic structure is constrained
% to be planar in the $xy$-plane. The [gm_planard] constraint function is
% used where the first 3 parameter determined the propagation vector,
% followed by the polar angles of the magnetic moments (here there is only
% 1 magnetic moment in the unit cell) which is fixed to 0. Finally the last
% 2 parameters corresponds to the polar angles of the normal to the
% spin-plane which is the $z$-axis ($\theta=0$, $\varphi=0$). The optimized
% magnetic structure is plotted.
%
% ```
% >>tri = sw_model('triAF',1)
% >>X1 = [0 0 0 0 0 0]
% >>X2 = [0 1/2 1/2 0 0 0]
% >>optRes = tri.optmagstr('func',@gm_planard,'xmin',X1,'xmax',X2)
% >>km = optRes.x(1:3)>>
% >>plot(tri)
% >>>swplot.zoom(1.5)
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
% `'func'`
% : Function that produces the spin orientations, propagation vector and
%   normal vector from the optimization parameters and has the following
%   argument list:
%   ```
%   [M, k, n] = @(x)func(M0, x)
%   ```
%  here `M` is matrix with dimensions of $[3\times n_{magExt}]$, `k` is the
%  propagation vector (row vector with 3 elements), `n` is the normal vector
%  of the spin rotation plane (row vector with 3 elements). The
%  default value is `@gm_spherical3d`. For planar magnetic structures
%  use `@gm_planar`.
% 
% `'xmin'`
% : Lower limit of the optimisation parameters.
% 
% `'xmax'`
% : Upper limit of the optimisation parameters.
% 
% `'x0'`
% : Starting value of the optimisation parameters. If empty
%   or undefined, then random values are used within the given limits.
% 
% `'boundary'`
% : Boundary conditions of the magnetic cell:
%   * `'free'`  Free, interactions between extedned unit cells are
%             omitted.
%   * `'per'`   Periodic, interactions between extended unit cells
%             are retained.
%
%   Default value is `{'per' 'per' 'per'}`.
% 
% `'epsilon'`
% : The smallest value of incommensurability that is tolerated
%   without warning. Default value is $10^{-5}$.
% 
% `'nRun'`
% : Number of runs. If random starting parameters are given, the
%   optimisation process will be rerun `nRun` times and the best
%   result (lowest ground state energy per spin) will be kept.
% 
% `'title'`
% : Gives a title string to the simulation that is saved in the
%   output.
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
% #### Limits on selected prameters
%
% Limits can be given on any input parameter of the constraint function by
% giving the name of the parameter. For parameter names see the help of the
% used constraint function. Limits per optimization parameter can be given
% in the following format: `optmagstr('ParName',[min max],...)`. For example
% to fix the `nTheta` value of [gm_planar] during the optimisation to zero
% use: `optmagstr(obj,'func',@gm_planar,'nTheta',[0 0])`.
%
% 
% #### Optimisation parameters
% 
% The optimization parameters are identical to the input options of the
% Matlab built-in optimizer [matlab.fminsearch].
%
% `'tolx'`
% : Minimum change of `x` when convergence reached, default
%     value is $10^{-4}$.
% 
% `'tolfun'`
% : Minimum change of the $R$ value when convergence reached,
%     default value is $10^{-5}$.
% 
% `'maxfunevals'`
% : Maximum number of function evaluations, default value
%     is $10^7$.
% 
% `'maxiter'`
% : Maximum number of iterations, default value is $10^4$.
% 
% ### Output Arguments
% 
% `optm`
% : Struct type variable with the following fields:
%   * `obj`       spinw object that contains the optimised magnetic structure.
%   * `x`         Optimised paramters in a row vector with $n_{par}$ number
%                 of elements.
%   * `fname`     Name of the contraint function.
%   * `xname`     Cell containing the name of the $x$ parameters with
%                   $n_{par}$ elements.
%   * `e`         Energy per spin in the optimised structure.
%   * `exitflag`  Exit flag of the optimisation code, see [matlab.fminsearch].
%   * `output`    Detailed output of the optimisation code, see [matlab.fminsearch].
%   * `param`     Input parameters, stored in a struct.
% 
% ### See Also
% 
% [spinw] \| [spinw.anneal] \| [gm_spherical3d] \| [gm_planar] \| [fminsearch]
%

if ~any(obj.atom.mag)
    error('spinw:optmagstr:NoMagAtom','There is no magnetic atom in the unit cell with S>0!');
end

% save the time of the beginning of the calculation
if nargout > 0
    optm.datestart  = datestr(now);
end

title0 = 'Optimised magnetic structure using simplex search';

inpForm.fname  = {'epsilon' 'func'           'boundary'          'xmin'   'xmax'  'x0'   };
inpForm.defval = {1e-5      @gm_spherical3d  {'per' 'per' 'per'} []       []      []     };
inpForm.size   = {[1 1]     [1 1]            [1 3]               [1 -1]   [1 -2]  [1 -3] };
inpForm.soft   = {0         0                0                   1        1       1      };

inpForm.fname  = [inpForm.fname  {'tolx' 'tolfun' 'maxfunevals' 'nRun' 'maxiter' 'title' 'tid'}];
inpForm.defval = [inpForm.defval {1e-4   1e-5     1e7           1      1e4       title0  -1   }];
inpForm.size   = [inpForm.size   {[1 1]  [1 1]    [1 1]         [1 1]  [1 1]     [1 -4]  [1 1]}];
inpForm.soft   = [inpForm.soft   {0      0        0             0      0         1       false}];

warnState = warning('off','sw_readparam:UnreadInput');
param = sw_readparam(inpForm, varargin{:});
pref = swpref;

% If magnetic structure hasn't been initialised
% just use nExt = [1 1 1]
if isempty(obj.mag_str.F)
    S = obj.unit_cell.S;
    nExt = [1 1 1];
else
    S = sqrt(sum(obj.magstr.S.^2,1));
    nExt = double(obj.magstr.N_ext);
end
nMagExt = length(S);

if param.tid == -1
    param.tid = pref.tid;
end


% determine the limits from the constraint function
xparam_warn_id = 'spinw:optmagstr:WrongLengthXParam';
xparam_warn_msg = ['Provided %s is the wrong length, expected length %i, got %i. ' ...
                   'This input parameter will be ignored.'];
if nargout(param.func) == 6
    [~, ~, ~, fname, pname, limit] = param.func(S,[]);
    % limits of the fitting parameters
    nPar = size(limit,2);
    % limits of the parameters
    if numel(param.xmin) ~= nPar
        if ~isempty(param.xmin)
            warning(xparam_warn_id, xparam_warn_msg, 'xmin', nPar, numel(param.xmin));
        end
        param.xmin = limit(1,:);
    end
    
    if numel(param.xmax) ~= nPar
        if ~isempty(param.xmax)
            warning(xparam_warn_id, xparam_warn_msg, 'xmax', nPar, numel(param.xmax));
        end
        param.xmax = limit(2,:);
    end
    
    % check if any parameter name is given as explicit input with limits
    inpForm.fname  = pname;
    inpForm.defval = repmat({[]},1,nPar);
    inpForm.size   = repmat({[1 2]},1,nPar);
    inpForm.soft   = repmat({1},1,nPar);
    % test input parameters
    fparam = sw_readparam(inpForm, varargin{:});
    warning(warnState);
    for ii = 1:nPar
        if ~isempty(fparam.(inpForm.fname{ii}))
            param.xmin(ii) = fparam.(inpForm.fname{ii})(1);
            param.xmax(ii) = fparam.(inpForm.fname{ii})(2);
        end
    end
    % restrict input parameter limits to the hardwired constraint function
    % limits
    param.xmin = max(param.xmin,limit(1,:));
    param.xmax = min(param.xmax,limit(2,:));
    
else
    if isempty(param.xmin) || isempty(param.xmax) || (numel(param.xmin) ~= numel(param.xmax))
        error('spinw:optmagtr:WrongInput','Missing limits on the x fitting parameters (use xmin and xmax options)!');
    end
    nPar  = numel(param.xmin);
    fname = '';
    pname = repmat({''},1,nPar);
end

% get magnetic couplings
[SS, SI] = obj.intmatrix;

% add dipolar interactions
SS.all = [SS.all SS.dip];

% Initial parameters are random if param.x0 is undefined/wrong size
xRand = (numel(param.x0)~= nPar);
if xRand && ~isempty(param.x0)
    warning(xparam_warn_id, xparam_warn_msg, 'x0', nPar, numel(param.x0));
end

if (~xRand) || (param.nRun<1)
    param.nRun = 1;
end

% Modify the interaction matrices according to the boundary conditions.
for ii = 1:3
    if strcmp('free',param.boundary{ii})
        SS.all(:,SS.all(ii,:)~=0) = [];
    end
end

dR    = [SS.all(1:3,:) zeros(3,nMagExt)];
atom1 = [SS.all(4,:)   1:nMagExt];
atom2 = [SS.all(5,:)   1:nMagExt];
JJ    = cat(3,reshape(SS.all(6:14,:),3,3,[]),SI.aniso);
g     = SI.g;

% B * g: field * g-tensor
Bg  = permute(mmat(SI.field,g)*obj.unit.muB,[2 3 1]);

minE = 0;
minX = zeros(1,nPar);

sw_timeit(0,1,param.tid,'Optimizing magnetic structure');

dx = param.xmax - param.xmin;

% Loop over nRun times
for ii = 1:param.nRun
    if xRand
        x0 = rand(1,nPar).*dx+param.xmin;
    else
        x0 = param.x0;
    end
    [X, E, exitflag, output] = sw_fminsearchbnd(@(x)efunc(x, S, dR, atom1, atom2, JJ, nExt, Bg, param.epsilon, param.func),x0,param.xmin,param.xmax,...
        optimset('TolX',param.tolx,'TolFun',param.tolfun,'MaxFunEvals',param.maxfunevals,'MaxIter',param.maxiter,'Display','off'));
    
    if E < minE
        minE = E;
        minX = X;
    end
    
    sw_timeit(ii/param.nRun*100,0,param.tid);
        
end

sw_timeit(100,2,param.tid);

[M, k, n] = param.func(S, minX);

%obj.mag_str.S = M;
%obj.mag_str.k = k;
%obj.mag_str.n = n;
obj.genmagstr('mode','helical','S',M,'k',k,'n',n);

spinw.validate(obj);

% Create output struct
if nargout > 0
    optm.obj      = copy(obj);
    optm.x        = minX;
    optm.e        = minE;
    optm.exitflag = exitflag;
    optm.output   = output;
    optm.param    = param;
    optm.fname    = fname;
    optm.xname    = pname;
    optm.dateend  = datestr(now);
    optm.title    = param.title;
end

end

function E = efunc(x, S, dR, atom1, atom2, JJ, nExt, Bg, epsilon, func)
% Energy function, cost function to optimize.

[M, k, n] = func(S, x);

M1 = M(:,atom1);
M2 = M(:,atom2);

kExt    = k.*nExt;
nMagExt = length(S);

% Incommenasurate structure in the extended unit cell.
incomm = any(abs(k-round(k))>epsilon);

% Rotate spins for incommensurate structures.
if incomm && any(any(dR))
    dRIdx = find(any(dR));
    for ii = 1:length(dRIdx)
        M2(:,dRIdx(ii)) = sw_rot(n, kExt*dR(:,dRIdx(ii))*2*pi, M2(:,dRIdx(ii)));
    end
end

Ml = repmat(permute(M1,[1 3 2]),[1 3 1]);
Mr = repmat(permute(M2,[3 1 2]),[3 1 1]);

E =  (sum(sum(sum(Ml.*JJ.*Mr))) - sum(Bg(:).*M(:)))/nMagExt;

end