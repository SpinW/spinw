function fitsp = fitpow(obj, data, varargin)
% fits experimental powder spin wave data
%
% ### Syntax
%
% `fitsp = fitpow(obj, data, Name, Value)`
%
% ### Description
% `fitsp = fitpow(obj, data, Name, Value)` takes a set of energy cuts of powder
% spin wave spectrum and fits this.
%
% ### Input Arguments
%
% `obj`
% : [spinw] object.
%
% `data`
% : input data cuts which may be an array of any of:
%   * a Horace 1D sqw or d1d object
%   * a struct with fields 'x', 'y', 'e', 'qmin', 'qmax'
%   * a struct with fields 'IX', 'qmin', 'qmax'
%   * a cell array {x y e qmin qmax}
%   * a cell array {IX_dataset_1d qmin qmax}
%   where x is a vector of energy transfer,
%         y is a vector of intensity,
%         e is a vector of uncertainties (errorbar)
%         qmin is a scalar denoting the lower Q value of the cut
%         qmax is a scalar denoting the upper Q value of the cut
%         IX / IX_dataset_1d is a Horace IX_dataset_1D object
%   Instead of the pair of scalars qmin, qmax, you specify a
%   2-vector: [qmin qmax] under the key 'qrange', e.g.
%   * a struct with fields 'x', 'y', 'e', 'qrange'
%   * a cell array {x y e [qmin qmax]}
%
% ### Name-Value Pair Arguments
%
% `'func'`
% : Function to change the Hamiltonian in `obj`, it needs to have the
%   following header:
%   ```
%   obj = @func(obj, x)
%   ```
%
% `'xmin'`
% : Lower limit of the optimisation parameters, optional.
%
% `'xmax'`
% : Upper limit of the optimisation parameters, optional.
%
% `'x0'`
% : Starting value of the optimisation parameters. If empty
%  or undefined, random values are used within the given limits.
%
% `'optimizer'`
% : String that determines the type of optimizer to use, possible values:
%   * `'lm3'`       (Default) Levenberg-Marquardt as implemented in Horace, see [ndbase.lm3],
%   * `'pso'`       Particle-swarm optimizer, see [ndbase.pso],
%   * `'simplex'`   Matlab built-in simplex optimizer, see [fminsearch](www.mathworks.ch/help/matlab/ref/fminsearch.html).
%
% `'nRun'`
% : Number of consecutive fitting runs, each result is saved in the
%   output `fitsp.x` and `fitsp.R` arrays. If the Hamiltonian given by the
%   random `x` parameters is incompatible with the ground state,
%   those `x` values will be omitted and new random `x` values will be
%   generated instead. Default value is 1.
%
% `'nMax'`
% : Maximum number of runs, including the ones that produce error
%   (due to incompatible ground state). Default value is 1000.
%
% `'nQ'`
% : Number of |Q| points to use for each cut to simulate cut thickness
%   Default is 10.
%
% `'hermit'`
% : Method for matrix diagonalization, for details see [spinw.spinwave].
%
% `'formfact'`
% : Whether to include the form factor in the calculation (Default: true)
%
% `'formfactfun'`
% : The function used to calculate the form factor (Default: sw_mff)
%
% `'conv_func'`
% : The convolution function used to broaden the data in |Q| and E.
%   Default: no convolution is applied.
%
% `'imagChk'`
% : Checks that the imaginary part of the spin wave dispersion is
%   smaller than the energy bin size.
%   If false, will not check
%   If 'penalize' will apply a penalty to iterations that yield imaginary modes
%   If true, will stop the fit if an iteration gives imaginary modes
%   Default is `penalize`.
%
% `'intensity_scale'`
% : The scale factor for the intensity. Use a negative number for a fixed scale factor,
%   positive for it to be fitted from the given initial value,
%   and 0 for it to be fitted with an initial value guessed from the data.
%   Default: -1  (fixed, no scale factor applied).
%
% `'background'`
% : The type of background to use. Either 'none', 'flat', 'linear', 'quadratic'
%
% `'background_par'`
% : The initial background parameters as a vector from highest order.
%   E.g. if using ('background', 'quadratic'), then ('background_par', [a, b, c])
%   where the background used will be a*x.^2 + b*x + c.
%   Default: []  (empty vector - initial parameters will be guessed from data).
%
% Parameters for monitoring the fitting:
%
% `'plot'`
% : If `true`, the data and fits are plotted as the fitting process progresses.
%   Default is `true`.
%
% `'verbose'`
% : Whether to print information as the fit progresses
%   Default is `true`.
%
% Optimizer options:
%
% `'TolX'`
% : Minimum change of` x` when convergence reached, default
%   value is $10^{-4}$.
%
% `'TolFun'`
% : Minimum change of the R value when convergence reached,
%   default value is $10^{-5}$.
%
% `'maxfunevals'`
% : Maximum number of function evaluations, default value is $10^3$.
%
% `'maxiter'`
% : Maximum number of iterations for the [ndbse.pso] optimizer.
%   Default value is 20.
%
% ### Output Arguments
%
% Output `fitsp` is struct type with the following fields:
% * `obj`   Copy of the input `obj`, with the best fitted
%           Hamiltonian parameters.
% * `x`     Final values of the fitted parameters, dimensions are
%           $[n_{run}\times n_{par}]$. The rows of `x` are sorted according
%           to increasing R values.
% * `redX2` Reduced $\chi^2_\eta$ value, goodness of the fit stored in a column
%           vector with $n_{run}$ number of elements, sorted in increasing
%           order. $\chi^2_\eta$ is defined as:
%
%   $\begin{align}
%                   \chi^2_\eta &= \frac{\chi^2}{\eta},\\
%                   \eta        &= n-m+1,
%   \end{align}$
%   where \\eta is the degree of freedom, $n$ number of
%   observations and $m$ is the number of fitted parameters.
%
% * `exitflag`  Exit flag of the `fminsearch` command.
% * `output`    Output of the `fminsearch` command.
%
% {{note As a rule of thumb when the variance of the measurement error is
% known a priori, \\chi$^2_\eta$\\gg 1 indicates a poor model fit. A
% \\chi$^2_\eta$\\gg 1 indicates that the fit has not fully captured the
% data (or that the error variance has been underestimated). In principle,
% a value of \\chi$^2_\eta$= 1 indicates that the extent of the match
% between observations and estimates is in accord with the error variance.
% A \\chi$^2_\eta$ < 1 indicates that the model is 'over-fitting' the data:
% either the model is improperly fitting noise, or the error variance has
% been overestimated.}}
%
% Any other option used by [spinw.spinwave] function are also accepted.
%
% ### See Also
%
% [spinw.spinwave] \| [spinw.matparser] \| [sw_egrid] \| [sw_neutron]
%

pref = swpref;
tid0 = pref.tid;

if nargin < 2
    error('spinw:fitpow:MissingData','No input data cuts given')
end

inpForm.fname  = {'xmin'   'xmax'  'x0'    'func' 'fibo' 'nRand' 'verbose'};
inpForm.defval = {[]       []      []      []     true   100     true};
inpForm.size   = {[1 -3]   [1 -4]  [1 -5]  [1 1]  [1 1]  [1 1]   [1 1]};
inpForm.soft   = {true     true    true    false  false  false   false};

inpForm.fname  = [inpForm.fname  {'formfact' 'formfactfun' 'conv_func' 'nQ'}];
inpForm.defval = [inpForm.defval {true       @sw_mff       @(s)s       10}];
inpForm.size   = [inpForm.size   {[1 -2]     [1 1]         [1 1]       [1 1]}];
inpForm.soft   = [inpForm.soft   {false      false         false       false}];

inpForm.fname  = [inpForm.fname  {'intensity_scale' 'background' 'background_par'}];
inpForm.defval = [inpForm.defval {-1                'linear'     []}];
inpForm.size   = [inpForm.size   {[1 1]             [1 -9]       [1 -10]}];
inpForm.soft   = [inpForm.soft   {false             false        true}];

inpForm.fname  = [inpForm.fname  {'tolx' 'tolfun' 'maxfunevals' 'nRun' 'plot'}];
inpForm.defval = [inpForm.defval {1e-4   1e-5     1e3           1      true}];
inpForm.size   = [inpForm.size   {[1 1]  [1 1]    [1 1]         [1 1]  [1 1]}];
inpForm.soft   = [inpForm.soft   {false  false    false         false  false}];

inpForm.fname  = [inpForm.fname  {'nMax' 'hermit' 'optimizer'}];
inpForm.defval = [inpForm.defval {1e3    true     'lm3'      }];
inpForm.size   = [inpForm.size   {[1 1]  [1 1]    [1 -7]     }];
inpForm.soft   = [inpForm.soft   {false  false    false      }];

inpForm.fname  = [inpForm.fname  {'maxiter' 'sw'  'optmem' 'tid' 'imagChk'}];
inpForm.defval = [inpForm.defval {20        1     0        tid0  'penalize'}];
inpForm.size   = [inpForm.size   {[1 1]  	[1 1] [1 1]    [1 1]  [1 -8]  }];
inpForm.soft   = [inpForm.soft   {false  	false false    false  false   }];

param = sw_readparam(inpForm, varargin{:});

% number of parameters (length of x)
nPar = max(max(length(param.xmin),length(param.xmax)),length(param.x0));
nRun = param.nRun;

% Parses the input data
[lindat, data] = parse_cuts(data);

% setup parameters for spinw.spinwave()
param0      = param;
%param0.plot = false;
param0.tid  = 0;

% Parses the imagChk parameter - we need to set it to true/false for spinw.spinwave()
if strncmp(param.imagChk, 'penalize', 1)
    param0.imagChk = true;
    param0.penalize_imag = true;
else
    param0.imagChk = logical(param.imagChk(1));
    param0.penalize_imag = false;
end

optfunname = sprintf('ndbase.%s', param.optimizer);
assert(~isempty(which(optfunname)), 'spinw:fitpow:invalidoptimizer', 'Optimizer %s not recognised', param.optimizer);
optimizer = str2func(optfunname);

%x     = zeros(nRun,nPar);
redX2 = zeros(nRun,1);

sw_timeit(0, 1, param.tid,'Fitting powder spin wave spectra');

idx = 1;

while idx <= nRun
    if ~isempty(param.x0)
        x0 = param.x0;
    elseif isempty(param.xmax) && isempty(param.xmin)
        x0 = rand(1,nPar);
    elseif isempty(param.xmax)
        x0 = param.xmin + rand(1,nPar);
    elseif isempty(param.xmin)
        x0 = param.xmax - rand(1,nPar);
    else
        x0 = rand(1,nPar).*(param.xmax-param.xmin)+param.xmin;
    end
    xmin = param.xmin;
    xmax = param.xmax;

    % Handle scale factor
    if param.intensity_scale < 0
        param0.fit_scale = false;
    else
        if param.intensity_scale == 0
            sf = guess_intensity(obj, data, param0);
        else
            sf = param.intensity_scale;
        end
        param0.fit_scale = true;
        x0 = [x0 sf];
        if ~isempty(xmin)
            xmin = [xmin 0];
        end
        if ~isempty(xmax)
            xmax = [xmax 10*sf];
        end
    end

    % Handle background
    param0.bkgfun = [];
    if ~strcmp(param.background, 'none')
        xb = param.background_par;
        if strcmp(param.background, 'flat')
            param0.bkgfun = @(x, p) ones(numel(x),1)*p;
            if isempty(xb)
                xb = min(lindat.y);
            end
            xnx = 0;
            xmx = max(lindat.y) / 2;
        elseif strcmp(param.background, 'linear')
            param0.bkgfun = @(x, p) x*p(1) + p(2);
            if isempty(xb)
                xb = guess_bkg(data, 'linear');
            end
            xnx = [-abs(xb(1))*2, 0];
            xmx = [abs(xb(1))*2, max(lindat.y)/2];
        elseif strcmp(param.background, 'quadratic')
            param0.bkgfun = @(x, p) x.^2*p(1) + x*p(2) + p(3);
            if isempty(xb)
                xb = guess_bkg(data, 'quadratic');
            end
            xnx = [-abs(xb(1:2))*2, 0];
            xmx = [abs(xb(1:2))*2, max(lindat.y)/2];
        else
            error('spinw:fitpow:unknownbackground', 'Background type %s not supported', param.background);
        end
        xb(isnan(xb)) = 0;
        x0 = [x0 xb]; %#ok<*AGROW>
        param0.nbkgpar = numel(xb);
        if ~isempty(xmin)
            xmin = [xmin xnx];
        end
        if ~isempty(xmax)
            xmax = [xmax xmx];
        end
    end

    % Only cache powder calcs for gradient methods
    if strncmp(param.optimizer, 'lm', 2)
        pow_fitfun('clear_cache');
    else
        pow_fitfun('no_cache');
    end
    [x(idx,:),~, output(idx)] = optimizer(lindat, @(x,p)pow_fitfun(obj, lindat, data, param.func, p, param0), x0, ...
        'lb', xmin, 'ub', xmax, 'TolX', param.tolx, 'TolFun', param.tolfun, 'MaxIter', param.maxiter);
    redX2(idx) = output.redX2;

    idx = idx + 1;
    sw_timeit(idx/nRun*100,0,param.tid);
end

sw_timeit(100,2,param.tid);

% Sort results
[redX2, sortIdx] = sort(redX2);
x = x(sortIdx,:);

% set the best fit to the spinw object
param.func(obj,x(1,1:nPar));

% Store all output in a struct variable.
fitsp.obj      = copy(obj);
fitsp.x        = x;
fitsp.redX2    = redX2;
%fitsp.exitflag = exitflag;
fitsp.output   = output;

end

function out = horace1d_tostruct(dat)
    assert(isa(dat, 'd1d') || isa(dat, 'sqw'), 'spinw:fitpow:invalidinput', 'Input must be a Horace 1D object, cell or struct');
    if isa(dat, 'sqw')
        assert(dat.dimensions == 1, 'spinw:fitpow:invalidinput', 'Input SQW object must be 1D');
        dd = dat.data;
    else
        dd = dat;
    end
    assert(strcmp(dat.ulabel{1}, '|Q|'), 'spinw:fitpow:invalidinput', 'Input Horace object is not a powder cut');
    out = struct('x', (dd.p{1}(1:end-1) + dd.p{1}(2:end))/2, 'y', dd.s, 'e', dd.e, 'qmin', dd.iint(1,1), 'qmax', dd.iint(2,1));
end

function out = verify_cut_struct(dat)
    if all(isfield(dat, {'x', 'y', 'e', 'qmin', 'qmax'}))
        out = dat;
    elseif all(isfield(dat, {'x', 'y', 'e', 'qrange'}))
        out = struct('x', dat.x, 'y', dat.y, 'e', dat.e, 'qmin', dat.qrange(1), 'qmax', dat.qrange(2));
    elseif all(isfield(dat, {'IX', 'qmin', 'qmax'}))
        out = struct('x', dat.IX.x, 'y', dat.IX.signal, 'e', dat.IX.error, 'qmin', dat.qmin, 'qmax', dat.qmax);
    elseif all(isfield(dat, {'IX', 'qrange'}))
        out = struct('x', dat.IX.x, 'y', dat.IX.signal, 'e', dat.IX.error, 'qmin', dat.qrange(1), 'qmax', dat.qrange(2));
    else
        error('spinw:fitpow:invalidinput', 'Input structure must have fields (x,y,e) or (IX) and (qmin,qmax) or (qrange)');
    end
end

function out = xyecell_tostruct(dat)
    errmsg = 'Cell input should consist of IX_dataset objects or numeric arrays';
    if numel(dat) > 3
        assert(all(cellfun(@isnumeric, dat)), 'spinw:fitpow:invalidinput', errmsg);
        xye = dat(1:3);
        remd = dat(4:end);
    elseif numel(dat) > 1
        assert(isa(dat{1}, 'IX_dataset_1d'), 'spinw:fitpow:invalidinput', errmsg);
        xye = {dat{1}.x, dat{1}.signal, dat{1}.error};
        remd = dat(2:end);
    end
    if numel(remainder) == 2
        out = struct('x', xye{1}, 'y', xye{2}, 'e', xye{3}, 'qmin', remd{1}, 'qmax', remd{2});
    else
        out = struct('x', xye{1}, 'y', xye{2}, 'e', xye{3}, 'qmin', remd{1}(1), 'qmax', remd{1}(2));
    end
end

function [lindat, outstruct] = parse_cuts(input)
    %outstruct = struct();
    if iscell(input)
        for ii = 1:numel(input)
            if isstruct(input{ii})
                outstruct(ii) = verify_cut_struct(input{ii});
            elseif iscell(input{ii})
                outstruct(ii) = xyecell_tostruct(input{ii});
            else
                outstruct(ii) = horace1d_tostruct(input{ii});
            end
        end
    elseif isstruct(input)
        for ii = 1:numel(input)
            outstruct(ii) = verify_cut_struct(input(ii));
        end
    else
        for ii = 1:numel(input)
            outstruct(ii) = horace1d_tostruct(input(ii));
        end
    end
    for ii = 1:numel(outstruct)
        xx = outstruct(ii).x(:).';
        df = diff(xx) ./ 2;
        outstruct(ii).evect = [xx(1)-df(1), xx(1:end-1) + df, xx(end) + df(end)];
    end
    lindat = linearise_data(outstruct);
end

function lindat = linearise_data(outstruct)
    % Linearise the data
    [x, y, e] = deal([], [], []);
    is_same_e = true;
    for ii = 1:numel(outstruct)
        x = [x; ii+([1:numel(outstruct(ii).y)].'/1e9)]; %#ok<NBRAK>
        y = [y; outstruct(ii).y(:)];
        e = [e; outstruct(ii).e(:)];
        if ii == 1
            evec = outstruct(ii).x;
        else
            if ~is_same_e && (numel(evec) ~= numel(outstruct(ii).x) ...
                             || sum(abs(evec(:) - outstruc(ii).x(:))) > 1e-5)
               is_same_e = false;
            end
        end
    end
    e(isnan(y)) = 1;
    y(isnan(y)) = 0;
    e(e==0) = 1;
    lindat = struct('x', x, 'y', y, 'e', e, 'is_same_e', is_same_e);
end

function out = is_same_data(d0, data)
    out = numel(d0) == numel(data) && numel(d0(1).x) == numel(data(1).x) && ...
        all(d0(1).x(:) == data(1).x(:));
end

function out = guess_bkg(data, bktype)
    xx = []; yy = [];
    for ii = 1:numel(data)
        xx = [xx; data(ii).x(:)];
        yy = [yy; data(ii).y(:)];
    end
    my = min(yy);
    idx = find(yy < (max(yy)-my)/20+my);
    [xx, isr] = sort(xx(idx));
    yy = yy(idx(isr));
    p2 = yy(end) / xx(end).^2;
    p1 = (yy(end) - yy(1)) / (xx(end) - xx(1));
    p0 = min(yy);
    if strcmp(bktype, 'linear')
        out = fminsearch(@(p) sum(abs(yy - (p(1)*xx + p(2)))), [p1 p0]);
    elseif strcmp(bktype, 'quadratic')
        out = fminsearch(@(p) sum(abs(yy - (p(1)*xx.^2 + p(2)*xx + p(3)))), [p2 p1 p0]);
    end
end

function out = guess_intscale(swobj, data, param)
    sfs = zeros(numel(data), 1);
    for ii = 1:numel(data)
        qq = linspace(data(ii).qmin, data(ii).qmax, param.nQ);
        spec = swobj.powspec(qq, 'Evect', data(ii).evect, 'nRand', 100, 'fibo', param.fibo, ...
            'hermit', param.hermit, 'formfact', param.formfact, 'formfactfun', param.formfactfun, ...
            'imagChk', param.imagChk);
        spec = param.conv_func(spec);
        yy = sum(spec.swConv, 2);
        sfs(ii) = max(data(ii).y) / max(yy);
    end
    out = mean(sfs);
end

function yCalc = pow_fitfun(obj, lindat, data, parfunc, x, param)
% calculates the powder spin wave spectrum to directly compare to data
%
% yCalc = SW_FITFUN(obj, data, swfunc, x, param)
%
% swobj         spinw object
% data          Structure contains the experimental data
% swfunc        Function to change the Hamiltonian in obj, of the form:
%                   obj = swfunc(obj,x);
% x             Actual parameter values.
% param         Parameters for the spinwave function.

persistent swp hf lines d0 useCache cache nCache iCache nEval;

if ischar(obj)
    if strcmp(obj, 'clear_cache')
        cache = [];
        nCache = [];
        iCache = 1;
        useCache = true;
    elseif strcmp(obj, 'no_cache')
        useCache = false;
    end
    nEval = 0;
    return;
end

if isempty(swp)
    swp = swpref;
end
fid0 = swp.fid;
swp.fid = 0;
cleanupObj = onCleanup(@()swpref.setpref('fid', fid0));

% Strip out background and intensity parameters before feeding the SpinW model
x_rem = x;
if param.verbose
    fprintf('Evaluation #%i with parameters:', nEval); fprintf('%f ', x_rem);
    nEval = nEval + 1;
end
sf = 1.0;
if param.intensity_scale < 0
    sf = abs(param.intensity_scale);
else
    if ~isempty(param.bkgfun)
        bkgpar = x_rem((end-param.nbkgpar+1):end);
        x_rem = x_rem(1:(end-param.nbkgpar));
    end
    if param.fit_scale
        sf = x_rem(end);
        x_rem = x_rem(1:end-1);
    end
end

% Modify the SpinW model with new parameters
parfunc(obj, x_rem);

% Other setup
fitstruc = data;
nD = numel(data);
nY = 0;

if param.plot
    if isempty(hf) || ~isvalid(hf) && (isempty(d0) || is_same_data(d0, data))
        hf = findobj('Tag', 'fitpow_plot');
    end
    if isempty(hf) || ~isvalid(hf)
        hf = figure('Tag', 'fitpow_plot');
        nR = ceil(nD / 4);
        for ii = 1:nD
            subplot(nR, ceil(nD / nR), ii); hold all;
            errorbar(data(ii).x, data(ii).y, data(ii).e, '.k');
            title(sprintf('%4.2f < |Q| < %4.2f', data(ii).qmin, data(ii).qmax));
            box on;
        end
        lines = [];
        d0 = data;
    end
    set(0, 'CurrentFigure', hf);
end

not_cached = true;
if ~isempty(cache) && useCache
    for ii = 1:numel(cache)
        if x_rem == cache{ii}{1}
            not_cached = false;
            fitstruc = cache{ii}{2};
            break
        end
    end
end
if not_cached
    if lindat.is_same_e
        nY = nD * numel(data(1).x);
        qq = zeros(nD * param.nQ);
        for ii = 1:nD
            i0 = (ii - 1) * param.nQ + 1;
            qq(i0:(i0 + param.nQ - 1)) = linspace(data(ii).qmin, data(ii).qmax, param.nQ);
        end
        spec = obj.powspec(qq, 'Evect', data(1).evect, 'nRand', param.nRand, 'fibo', param.fibo, ...
            'hermit', param.hermit, 'formfact', param.formfact, 'formfactfun', param.formfactfun, ...
            'imagChk', param.imagChk);
        spec = param.conv_func(spec);
        for ii = 1:nD
            i0 = (ii - 1) * param.nQ + 1;
            fitstruc(ii).y = sum(spec.swConv(:, i0:(i0 + param.nQ - 1)), 2);
        end
    else
        for ii = 1:nD
            qq = linspace(data(ii).qmin, data(ii).qmax, param.nQ);
            spec = obj.powspec(qq, 'Evect', data(ii).evect, 'nRand', param.nRand, 'fibo', param.fibo, ...
                'hermit', param.hermit, 'formfact', param.formfact, 'formfactfun', param.formfactfun, ...
                'imagChk', param.imagChk);
            spec = param.conv_func(spec);
            fitstruc(ii).y = sum(spec.swConv, 2);
            nY = nY + numel(fitstruc(ii).y);
        end
    end
    if useCache
        if isempty(nCache)  % Determines cache size from free memory up to max 1GB
            nCache = min([numel(x_rem)*5, ceil((min([sw_freemem, 1e9]) / 4) / (numel(fitstruc) * numel(fitstruc(1).y) * 3))]);
        end
        cache{iCache} = {x_rem, fitstruc};
        iCache = iCache + 1;
        if iCache > nCache
            iCache = 1;
        end
    end
end

yCalc = zeros(nY, 1);
i0 = 1;
for ii = 1:numel(data)
    %if param.fit_separate
    %    sf = 1.0;
    %    ssf = fminsearch(@(x)sum(abs(data(ii).y(:) - x*fitstruc(ii).y(:)), 'omitnan'), max(data(ii).y) / max(fitstruc(ii).y));
    %    if ssf < 1e5
    %        fitstruc(ii).y = ssf * fitstruc(ii).y;
    %    end
    %end
    if sf ~= 1.0
        fitstruc(ii).y = fitstruc(ii).y * sf;
    end
    i1 = i0 + numel(fitstruc(ii).y) - 1;
    if ~isempty(param.bkgfun)
        yCalc(i0:i1) = fitstruc(ii).y + param.bkgfun(data(ii).x, bkgpar);
    else
        yCalc(i0:i1) = fitstruc(ii).y;
    end
    i0 = i1 + 1;
end

%if param.fit_together
%    sf = fminsearch(@(x)sum(abs(lindat.y(:) - x*yCalc(:)), 'omitnan'), max(l0.y) / max(yCalc));
%end

yCalc(isnan(yCalc)) = 0;
if ~param.imagChk
    yCalc = abs(yCalc);
end

if param.plot
    nR = ceil(nD / 4);
    i0 = 1;
    for ii = 1:numel(data)
        subplot(nR, ceil(nD / nR), ii);
        try %#ok<TRYNC>
            delete(lines(ii));
        end
        i1 = i0 + numel(fitstruc(ii).y) - 1;
        lines(ii) = plot(fitstruc(ii).x, yCalc(i0:i1), '-r');
        i0 = i1 + 1;
    end
    drawnow;
end

if param.verbose
    fprintf('. chi2 = %d\n', sum((yCalc - lindat.y).^2 ./ lindat.e.^2) / numel(yCalc));
end

end
