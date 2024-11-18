function [pOpt,fVal,stat] = simplex(dat,func,p0,varargin)
% simplex optimisation
%
% ### Syntax
%
% `[pOpt,fVal,stat] = ndbase.simplex([],func,p0,Name,Value)`
%
% `[pOpt,fVal,stat] = ndbase.simplex(dat,func,p0,Name,Value)`
%
% ### Description
% 
% `[pOpt,fVal,stat] = ndbase.simplex([],func,p0,Name,Value)` finds a
% minimum of a function of several parameters using the constrained simplex
% optimization algorithm by calling the unconstrained Matlab built-in
% algorithm [matlab.fminsearch].
%
% The function has two different modes, depending on the first input
% argument. If `dat` is empty, `simplex` minimizes the cost function func,
% that has the following header:
% ```
% R2 = func(p)
% ```
%
% If `dat` is a struct, `simplex` optimizes the model defined by `func` via
% least squares to the data stored in `dat`. In this case `func` has the
% following header:
% ```
% yModel = func(x,p)
% ```
%
% And the least square deviation is defined by:
%
% $R^2 = \sum \frac{(y_{dat}-y_{fit})^2}{\sigma_{dat}^2}$
%
%  {{note If options is supplied, then TolX will apply to the transformed
%  variables. All other [matlab.fminsearch] parameters should be unaffected.
%
%  Variables which are constrained by both a lower and an upper
%  bound will use a $\sin(x)$ transformation. Those constrained by
%  only a lower or an upper bound will use a quadratic
%  transformation, and unconstrained variables will be left alone.
%
%  Variables may be fixed by setting their respective bounds equal.
%  In this case, the problem will be reduced in size for [matlab.fminsearch].
%
%  The bounds are inclusive inequalities, which admit the
%  boundary values themselves, but will not permit any function
%  evaluations outside the bounds. These constraints are strictly
%  followed.
%
%  If your problem has an exclusive (strict) constraint which will
%  not admit evaluation at the bound itself, then you must provide
%  a slightly offset bound. An example of this is a function which
%  contains the log of one of its parameters. If you constrain the
%  variable to have a lower bound of zero, then `simplex` may
%  try to evaluate the function exactly at zero.}}
%
% ### Examples
%
% Example usage on the rosen function.
%
% ```
% >>rosen = @(x) (1-x(1)).^2 + 105*(x(2)-x(1).^2).^2
% ```
%
% Unconstrained optimisation:
%
% ```
% >>>rosen = @(x) (1-x(1)).^2 + 105*(x(2)-x(1).^2).^2
% >>fminsearch(rosen,[3 3])>>
% ```
%
% Constrained optimisation:
%
% ```
% >>>rosen = @(x) (1-x(1)).^2 + 105*(x(2)-x(1).^2).^2
% >>ndbase.simplex([],rosen,[3 3],'lb',[2 2],'ub',[])>>
% ```
%
% ### Input Arguments
%
% `dat`
% : Either empty or contains data to be fitted stored in a structure with
%   fields:
%   * `dat.x`   vector of $N$ independent variables,
%   * `dat.y`   vector of $N$ data values to be fitted,
%   * `dat.e`   vector of $N$ standard deviation (positive numbers)
%               used to weight the fit. If zero or missing
%               `1/dat.y^2` will be assigned to each point.
%
% `func`
% : Function handle with one of the following definition:
%   * `R = func(p)`         if `dat` is empty,
%   * `y  = func(x,p)`      if `dat` is a struct.
%   Here `x` is a vector of $N$ independent variables, `p` are the
%   $M$ parameters to be optimized and `y` is the simulated model.
%   If `resid_handle` argument is false (default) then the function returns 
%   a scalar (the cost funciton to minimise e.g. chi-squared). If 
%   `resid_handle` is true then the function returns a vector of residuals
%   (not the residuals squared).
%
% `p0`
% : vector of initial parameter guesses - starting point for the
%  optimisation.
%
% ### Name-Value Pair Arguments
%
% `'lb'`
% : Vector with $N$ elements, lower boundary of the parameters. Default
%   value is -inf.
%
% `'ub'`
% : Vector with $N$ elements, upper boundary of the parameters. Default
%   value is inf.
%
% `'MaxIter'`
% : Maximum number of iterations, default value is $100M$.
%
% `'MaxFunEvals'`
% : Maximum number of function evaluations, default value is
%   $1000M$. NOT IMPLEMENTED!
%
% `'TolX'`
% : Convergence tolerance for parameters, default value is $10^{-3}$.
%
% `'Display'`
% : Level of information to print onto the Command Line during
%   optimisation. Possible values are `'off'`, `'notify'` and `'iter'`.
%   Default value is `'off'`.
%
% `'TolFun'`
% : Convergence tolerance on the `R2` value (return value of `func` or
%   the weighted least square deviation from data). Default value is
%   $10^{-3}$.
%
% `'resid_handle'`
% : Boolean scalar - if true and `dat` is empty then 'func' fucntion handle 
%   returns array of residuals, if false (default) then function handle 
%   returns a scalar cost function.
%
% `'vary'`
% : Boolean vector with $M$ elements (one per parameter), if an element is 
%   false, the corresponding parameter will be fixed (not optimised). 
%   Default is true for all parameters.
%
% ### See Also
%
% [ndbase.lm] \| [ndbase.pso]
%
% Original author: John D'Errico
% E-mail: woodchips@rochester.rr.com
% Release: 4
% Release date: 7/23/06

    if nargin == 0
        swhelp ndbase.simplex
        return
    end
    
    % number of parameters
    Np = numel(p0);

    % not implemented yet: 'MaxFunEvals'
    inpForm.fname  = {'Display' 'TolFun' 'TolX' 'MaxIter' 'lb'      'ub'    };
    inpForm.defval = {'off'     1e-3     1e-3   100*Np    []        []      };
    inpForm.size   = {[1 -1]    [1 1]    [1 1]  [1 1]     [-5 -2]   [-3 -4] };
    inpForm.soft   = {false     false    false  false     true      true    };

    inpForm.fname  = [inpForm.fname  {'resid_handle', 'vary'}];
    inpForm.defval = [inpForm.defval {false,          true(1, Np)}];
    inpForm.size   = [inpForm.size   {[1 1],          [1, Np]}];
    
    param = sw_readparam(inpForm, varargin{:});
  
    cost_func_wrap = ndbase.cost_function_wrapper(func, p0, "data", dat, ...
                                              'lb', param.lb, 'ub', param.ub, ...
                                              'resid_handle', param.resid_handle, ...
                                              'ifix', find(~param.vary));
    
    % transform starting values into their unconstrained surrogates.
    p0_free = cost_func_wrap.get_free_parameters(p0);
    
    if isempty(p0_free)
        % All parameters fixed, evaluate cost at initial guess
        % don't use p0 as could contain fixed params outside bounds
        pOpt = cost_func_wrap.get_bound_parameters(p0_free);
        fVal = cost_func_wrap.eval_cost_function(p0_free);
        fit_stat.message = 'Parameters are fixed, no optimisation';
        fit_stat.iterations = 0;
        fit_stat.funcCount  = 1;
        exitFlag = 0;
    else
        % now we can call fminsearch, but with our own free parameter
        [p_free, fVal, exitFlag, fit_stat] = fminsearch(@cost_func_wrap.eval_cost_function, p0_free, param);
        % undo the variable transformations into the original space
        pOpt = cost_func_wrap.get_bound_parameters(p_free);
    end

    % setup output struct
    stat.sigP       = [];
    stat.Rsq        = [];
    stat.sigY       = [];
    stat.corrP      = [];
    stat.cvgHst     = [];
    stat.algorithm  = 'Nelder-Mead simplex direct search';
    stat.func = cost_func_wrap.cost_func;
    stat.param = param;
    stat.param.Np = cost_func_wrap.get_num_free_parameters();
    stat.msg        = fit_stat.message;
    stat.iterations = fit_stat.iterations;
    stat.funcCount  = fit_stat.funcCount;
    stat.exitFlag   = exitFlag;
    stat.p = pOpt;
    if isempty(dat)
        stat.redX2 = fVal;
    else
        % divide R2 with the statistical degrees of freedom
        stat.redX2   = fVal/(numel(dat.x)-stat.param.Np+1);
    end
    
end

