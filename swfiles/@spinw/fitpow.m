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
%   * `'pso'`       Particle-swarm optimizer, see [ndbase.pso],
%                   default.
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
% `'hermit'`
% : Method for matrix diagonalization, for details see [spinw.spinwave].
% 
% `'epsilon'`
% : Small number that controls wether the magnetic structure is
%   incommensurate or commensurate, default value is $10^{-5}$.
% 
% `'imagChk'`
% : Checks that the imaginary part of the spin wave dispersion is
%   smaller than the energy bin size. 
%   If false, will not check
%   If 'penalize' will apply a penalty to iterations that yield imaginary modes
%   If true, will stop the fit if an iteration gives imaginary modes
%   Default is `penalize`.
% 
% Parameters for visualizing the fit results:
% 
% `'plot'`
% : If `true`, the measured dispersion is plotted together with the
%   fit. Default is `true`.
% 
% `'iFact'`
% : Factor of the plotted simulated spin wave intensity (red
%   ellipsoids).
% 
% `'lShift'`
% : Vertical shift of the `Q` point labels on the plot.
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
% `'MaxFunEvals'`
% : Maximum number of function evaluations, default value is
%   $10^7$.
% 
% `'MaxIter'`
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

inpForm.fname  = {'epsilon' 'xmin'   'xmax'  'x0'    'func' 'fibo' 'nRand'};
inpForm.defval = {1e-5      []       []      []      []     true   100};
inpForm.size   = {[1 1]     [1 -3]   [1 -4]  [1 -5]  [1 1]  [1 1]  [1 1]};
inpForm.soft   = {true      true     true    true    false  false  false};

inpForm.fname  = [inpForm.fname  {'formfact' 'formfactfun' 'conv_func' 'nQ'  'intensity_scale'}];
inpForm.defval = [inpForm.defval {false      @sw_mff       @(s)s       10    0}];
inpForm.size   = [inpForm.size   {[1 -2]     [1 1]         [1 1]       [1 1] [1 -9]}];
inpForm.soft   = [inpForm.soft   {false      false         false       false true}];

inpForm.fname  = [inpForm.fname  {'tolx' 'tolfun' 'maxfunevals' 'nRun' 'plot'}];
inpForm.defval = [inpForm.defval {1e-4   1e-5     1e3           1      true}];
inpForm.size   = [inpForm.size   {[1 1]  [1 1]    [1 1]         [1 1]  [1 1]}];
inpForm.soft   = [inpForm.soft   {false  false    false         false  false}];

inpForm.fname  = [inpForm.fname  {'nMax' 'hermit' 'iFact' 'lShift' 'optimizer'}];
inpForm.defval = [inpForm.defval {1e3    true     1/30     0       'pso'      }];
inpForm.size   = [inpForm.size   {[1 1]  [1 1]    [1 1]    [1 1]   [1 -7]     }];
inpForm.soft   = [inpForm.soft   {false  false    false    false   false      }];

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

% Handle scale factor
param0.fit_separate = false;
param0.fit_together = false;
param0.scale_factor = 1.0;
if ischar(param.intensity_scale) || isstring(param.intensity_scale)
    if strcmp(param.intensity_scale, 'fit_separate')
        param0.fit_separate = true;
    elseif strcmp(param.intensity_scale, 'fit_together')
        param0.fit_together = true;
    end
elseif isfloat(param.intensity_scale)
    param0.scale_factor = param.intensity_scale;
end

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

x     = zeros(nRun,nPar);
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

    [x(idx,:),~, output(idx)] = optimizer(lindat, @(x,p)pow_fitfun(obj, lindat, data, param.func, p, param0), x0, ...
        'lb', param.xmin, 'ub', param.xmax, ...
        'TolX', param.tolx, 'TolFun', param.tolfun, 'MaxIter', param.maxiter);
    redX2(idx) = output.redX2;
    
    idx = idx + 1;
    sw_timeit(idx/nRun*100,0,param.tid);
end

sw_timeit(100,2,param.tid);

% Sort results
[redX2, sortIdx] = sort(redX2);
x = x(sortIdx,:);

% set the best fit to the spinw object
param.func(obj,x(1,:));

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
    out = struct('x', dd.p{1}, 'y', dd.s, 'e', dd.e, 'qmin', dd.iint(1,1), 'qmax', dd.iint(2,1));
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
        x = [x; ii+([1:numel(outstruct(ii).y)]/1e9)];
        y = [y; outstruct(ii).y];
        e = [e; outstruct(ii).e];
        if ii == 1
            evec = outstruct(ii).x;
        else
            if ~is_same_e && (numel(evec) ~= numel(outstruct(ii).x) ... 
                             || sum(abs(evec(:) - outstruc(ii).x(:))) > 1e-5)
               is_same_e = false;
            end
        end
    end
    lindat = struct('x', x, 'y', y, 'e', e, 'is_same_e', is_same_e);
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

persistent swp hf lines;
if isempty(swp)
    swp = swpref;
end
fid0 = swp.fid;
swp.fid = 0;
cleanupObj = onCleanup(@()swpref.setpref('fid', fid0));

parfunc(obj,x);

if param.plot
    if isempty(hf) || ~isvalid(hf)
        hf = findobj('Tag', 'fitpow_plot');
    end
    if isempty(hf) || ~isvalid(hf)
        hf = figure('Tag', 'fitpow_plot');
        for ii = 1:numel(data)
            subplot(1, numel(data), ii); hold all;
            errorbar(data(ii).x, data(ii).y, data(ii).e, '.k');
        end
        lines = [];
    end
    set(0, 'CurrentFigure', hf);
end

fitstruc = data;
sf = param.scale_factor;
nD = numel(data);
nY = 0;

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

yCalc = zeros(nY, 1);
i0 = 1;
for ii = 1:numel(data)
    if param.fit_separate
        sf = 1.0;
        ssf = fminsearch(@(x)sum(abs(data(ii).y(:) - x*fitstruc(ii).y(:)), 'omitnan'), max(data(ii).y) / max(fitstruc(ii).y));
        if ssf < 1e5
            fitstruc(ii).y = ssf * fitstruc(ii).y;
        end
    end
    i1 = i0 + numel(fitstruc(ii).y) - 1;
    yCalc(i0:i1) = fitstruc(ii).y;
    i0 = i1 + 1;
end

if param.fit_together
    sf = fminsearch(@(x)sum(abs(lindat.y(:) - x*yCalc(:)), 'omitnan'), max(l0.y) / max(yCalc));
end

if sf ~= 1.0
    yCalc = sf * yCalc;
end

if param.plot
    for ii = 1:numel(data)
        subplot(1, numel(data), ii);
        try %#ok<TRYNC>
            delete(lines(ii));
        end
        lines(ii) = plot(fitstruc(ii).x, sf * fitstruc(ii).y, '-r');
    end
    drawnow;
end

end
