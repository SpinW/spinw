function [pOpt, fVal, stat] = lm4(dat, func, p0, varargin)
% optimization of parameters using the Levenberg Marquardt method
%
% ### Syntax
%
% `[pOpt,fVal,stat] = ndbase.simplex([],func,p0,Name,Value)`
%
% `[pOpt,fVal,stat] = ndbase.simplex(dat,func,p0,Name,Value)`
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
% : Vector with $M$ elements, lower boundary of the parameters. Default
%   value is -inf (i.e. unbounded).
%
% `'ub'`
% : Vector with $M$ elements, upper boundary of the parameters. Default
%   value is inf (i.e. unbounded).
%
% `'resid_handle'`
% : Boolean scalar - if true and `dat` is empty then 'func' fucntion handle 
%   returns array of residuals, if false (default) then function handle 
%   returns a scalar cost function.
%
% `'diff_step'`
% : Vector with $M$ or 1 element, defines the fractional increment of
%   a parameter when calculating the Jacobians using the forward finite 
%   difference (default is 1e-7).
%
% `'MaxIter'`
% : Maximum number of iterations, default value is $100M$.
%
% `'gtol'`
% : Convergence tolerance on gradient vector of cost-function wrt change
%   in parameters (default 1e-8).
%
% `'ftol'`
% : Convergence tolerance on change in cost function (default 1e-8).
%
% `'ptol'`
% : Convergence tolerance on relative length of parameter step vector wrt
%   parameter vector (default 1e-8).
%
% `'lambda0'`
% : Initial Levenberg–Marquardt damping parameter which determines step 
%   size and angle of trajectory between gradient-descent (initially) and 
%   Gauss-Newton (where cost-function surface becomes increasingly 
%   parabolic). Initially `large` lambda values correponds to smaller steps
%   along gradient descent trajectory.  Decreasing `'lambda0'` may lead to 
%   faster convergence if the solution is near the minimum but will be less
%   reliable further from the minimum where the second-order expansion is
%   less valid. Default value is 1e-2.
%
% `'nu_up'`
% : Factor by which to scale the damping parameter (lambda) if parameter 
%   step increases the cost value. Typically > 1 (default is 5) - ie. 
%   increasing lambda to produce smaller steps such that the 
%   first-order approx. is more valid.
%
% `'nu_dn'`
% : Factor by which to scale the damping parameter (lambda) if parameter 
%   step decreases the cost value. Typically < 1 (default is 0.3) - i.e.
%   decreasing lambda to speed-up convergence by taking larger steps as 
%   the cost-function surface becomes increasingly parabolic.
%
% `'vary'`
% : Boolean vector with $M$ elements (one per parameter), if an element is 
%   false, the corresponding parameter will be fixed (not optimised). 
%   Default is true for all parameters.
%
% ### Output
%
% `pOpt`
% : Value of the $M$ optimal parameters.
%
% `fVal`
% : Value of the cost function calculated at the optimal parameters
%
% `stat`
% : Structure storing the detailed output of the calculation with
%  the following fields:
%   p         Least-squares optimal estimate of the parameter values.
%   redX2     Reduced Chi squared statistic, its value
%             should be close to 1. If the value is larger, the
%             model is not a good description of the data. If the
%             value is smaller, the model is overparameterized
%             and fitting the statistical error of the data.
%   sigP      Asymptotic standard error of the parameters.
%   sigY      Asymptotic standard error of the curve-fit (to implement).
%   corrP     Correlation matrix of the parameters (to implement).
%   Rsq       R-squared cofficient of multiple determination (to implement).
%   cvgHst    Convergence history (to implement)..
%   exitFlag  The reason, why the code stopped:
%               0       maximum number of iterations reached
%               1       convergence step-size (see `ptol`)
%               2       convergence in cost (see `ftol`)
%               3       convergence in gradient (see `gtol`)
%   msg       String, one of the above messages.
%   nIter     The number of iterations executed during the fit.
%   param     Input parameters passed to ndbase.lm4
%   algorithm name of algorithm (Levenberg Marquardt);
%   func      Function handle corresponding to cost-function
%
% See also ndbase.simplex.

nparams = numel(p0);

inpForm.fname  = {'diff_step' 'lb'            'ub'           'MaxIter' };
inpForm.defval = {1e-7        -inf(1,nparams) inf(1,nparams) 100*nparams};
inpForm.size   = {[1 -1]      [1 nparams]     [1 nparams]    [1 1]};

inpForm.fname  = [inpForm.fname  {'gtol' 'ftol' 'ptol' 'lambda0'}];
inpForm.defval = [inpForm.defval {1e-8    1e-8   1e-8  1e-2}];
inpForm.size   = [inpForm.size   {[1 1]   [1 1]  [1 1] [1 1]}];

inpForm.fname  = [inpForm.fname  {'nu_up', 'nu_dn',  'resid_handle'}];
inpForm.defval = [inpForm.defval {5        0.3,      false}];
inpForm.size   = [inpForm.size   {[1 1]    [1 1],    [1 1]}];

inpForm.fname  = [inpForm.fname  {'vary'}];
inpForm.defval = [inpForm.defval {true(1, nparams)}];
inpForm.size   = [inpForm.size   {[1 nparams]}];

param = sw_readparam(inpForm, varargin{:});

cost_func_wrap = ndbase.cost_function_wrapper(func, p0, "data", dat, ...
                                              'lb', param.lb, 'ub', param.ub, ...
                                              'resid_handle', param.resid_handle, ...
                                              'ifix', find(~param.vary));
% transform starting values into their unconstrained surrogates.
p0_free = cost_func_wrap.get_free_parameters(p0);

if isempty(dat) && ~param.resid_handle
    % minimising scalar - empty resid
    eval_cost_func = @eval_cost_scalar;
    calc_hessian_and_jacobian = @calc_hessian_and_jacobian_scalar;
    diff_step = param.diff_step; % always interpreted as fractional step even if 1 parameter
else
    % minimising sum square residuals
    eval_cost_func = @eval_cost_resids;
    calc_hessian_and_jacobian = @calc_hessian_and_jacobian_resid;
    % get absolute diff_step for each parameter
    diff_step = abs(p0_free).*param.diff_step;
    min_step = sqrt(eps);
    diff_step(abs(diff_step) < min_step) = min_step;
end

% eval at starting guess
[cost_val, resids] = eval_cost_func(cost_func_wrap, p0_free(:));
if isempty(resids)
    ndof = 1;  % minimising scalar function
else
    ndof = numel(resids) - numel(p0_free) + 1;
end
if isempty(p0_free)
    % Re-calc bound params as p0 as could have fixed params outside bounds
    pOpt = cost_func_wrap.get_bound_parameters(p0_free);
    fVal = cost_val/ndof;
    message = 'Parameters are fixed, no optimisation';
    perr = zeros(size(pOpt));
    niter = 0;
    exit_flag = 0;
else
    p = p0_free(:);
    lambda = param.lambda0;
    [hess, jac] = calc_hessian_and_jacobian_scalar(cost_func_wrap, p, diff_step, cost_val, resids);
    exit_flag = 0;
    for niter = 1:param.MaxIter
        dp = calc_parameter_step(hess, jac, lambda);
        if norm(dp) < param.ptol*(param.ptol + norm(p))
            message = "step size below tolerance ptol";
            exit_flag = 1;
            break
        end
        new_p = p - dp;
        new_cost_val = cost_func_wrap.eval_cost_function(new_p);
        dcost = new_cost_val - cost_val;
        if dcost < 0
            lambda = param.nu_dn*lambda; % decrease lambda
            p = new_p;
            [cost_val, resids] = eval_cost_func(cost_func_wrap, p);
            [hess, jac] = calc_hessian_and_jacobian(cost_func_wrap, p, diff_step, cost_val, resids);
            if abs(dcost) < param.ftol*ndof
                message = "change in reduced chi-sq below tolerance ftol";
                exit_flag = 2;
                break
            elseif norm(jac) < param.gtol
                message = "change in gradient vector below tolerance gtol";
                exit_flag = 3; 
                break
            end
        else
            % increase lambda so take path closer to steepest descent
            % don't recalc hess or jac
            lambda = param.nu_up*lambda;
        end
    end
    % collect output
    pOpt = cost_func_wrap.get_bound_parameters(p);
    perr = zeros(size(pOpt));
    fVal = cost_val / ndof;
    if exit_flag > 0
        % converged on solution - calculate errors
        cov = pinv(hess) * 2.0 * fVal;
        perr(cost_func_wrap.ifree) = sqrt(diag(cov));
    else
        message = "Failed to converge in MaxIter";
    end
end

stat.Rsq        = [];
stat.sigY       = [];
stat.corrP      = [];
stat.cvgHst     = [];
stat.algorithm  = 'Levenberg Marquardt';
stat.func = cost_func_wrap.cost_func;
stat.param = param;
stat.param.Np = cost_func_wrap.get_num_free_parameters();
stat.msg        = message;
stat.iterations = niter;
stat.exitFlag   = exit_flag;
stat.p = pOpt;
stat.redX2 = fVal;
stat.sigP = perr;

end

function dp = calc_parameter_step(hess, jac, lambda)
    damped_hess = hess + lambda*diag(diag(hess)); % LM
    if rcond(damped_hess) < eps
        dp = pinv(damped_hess)*jac; % nearly singular
    else
        try
            dp = damped_hess\jac;
        catch
            dp = pinv(damped_hess)*jac;
        end
    end
end

function [hess, jac] = calc_hessian_and_jacobian_scalar(cost_func_wrap, p, diff_step, cost_val, ~)
    % last input ignored (empty resids vector for consistent API)
    [hess, extra_info] = ndbase.estimate_hessian(@cost_func_wrap.eval_cost_function, p,  'cost_val', cost_val, 'step', diff_step);
    jac = extra_info.jacobian(:);
end

function [cost_val, resids] = eval_cost_scalar(cost_func_wrap, p)
    cost_val = cost_func_wrap.eval_cost_function(p);
    resids = [];
end

function [hess, jac] = calc_hessian_and_jacobian_resid(cost_func_wrap, p, diff_step, ~, resids)
    % evaluate jacobian of residuals using finite difference
    jac_resids = ones([numel(resids), numel(p)]);
    for ipar = 1:numel(p)
        p(ipar) = p(ipar) + diff_step(ipar);
        resids_one_step = cost_func_wrap.eval_resid(p);
        jac_resids(:,ipar) = (resids_one_step - resids)/diff_step(ipar);
        p(ipar) = p(ipar) - diff_step(ipar);
    end
    % eval jacobian of cost function 
    jac = 2*(jac_resids')*resids;
    % approx. hessian of cost fucntion
    hess = 2*(jac_resids')*jac_resids;
end

function [cost_val, resids] = eval_cost_resids(cost_func_wrap, p)
    resids = cost_func_wrap.eval_resid(p);
    cost_val = sum(resids(:).^2);
end

