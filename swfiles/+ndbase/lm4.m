function [pOpt, fVal, stat] = lm4(dat, func, p0, varargin)
% optimization of parameters using the Levenberg Marquardt method
%
% [pOpt,fVal,stat] = NDBASE.LM([],func,p0,'Option1','Value1',...)
%
% Levenberg Marquardt curve-fitting function, minimizes the return value of
% a given function.
%
% Input:
%
% func      Function handle with the following definition:
%               R = func(p)
%           where p are the M parameters to be optimized.
% p0        Vector of M initial parameters.
%
% Options:
%
% Options can be given using the modified output of optimset() or as option
% name string option value pairs.
%
% dp        Vector with N or 1 element, defines the fractional increment of
%           'p' when calculating the Jacobian matrix dFunc(x,p)/dp:
%               dp(j)>0     central differences calculated,
%               dp(j)<0     one sided 'backwards' differences calculated,
%               dp(j)=0     sets corresponding partials to zero, i.e. holds
%                           p(j) fixed.
%           Default value if 1e-3.
% vary      Vector with N elements, if an element is false, the
%           corresponding parameter will be fixed. Default value is
%           false(1,N).
% lb        Vector with N elements, lower boundary of the parameters.
%           Default value is -inf.
% ub        Vector with N elements, upper boundary of the parameters.
%           Default value is inf.
% MaxIter   Maximum number of iterations, default value is 100*M.
% TolX      Convergence tolerance for parameters, defines the maximum of
%           the relative chande of any parameter value. Default value is
%           1e-3.
% eps1      Convergence tolerance for gradient, default value is 1e-3.
% eps2      Convergence tolerance for reduced Chi-square, default value is
%           1e-2.
% eps3      Determines acceptance of a L-M step, default value is 0.1.
% lambda0   Initial value of L-M paramter, default value is 1e-2.
% nu0       Value that determines the speed of convergence. Default value
%           is 10. It should be larger than 1.
%
% Output:
%
% pOpt      Value of the M optimal parameters.
% fVal      Value of the model function calculated with the optimal
%           parameters at the N independent values of x.
%
% stat      Structure, storing the detailed output of the calculation with
%           the following fields:
%               p       Least-squares optimal estimate of the parameter
%                       values.
%               redX2   Reduced Chi squared error criteria, its value
%                       should be close to 1. If the value is larger, the
%                       model is not a good description of the data. If the
%                       value is smaller, the model is overparameterized
%                       and fitting the statistical error of the data.
%               sigP    Asymptotic standard error of the parameters.
%               sigY    Asymptotic standard error of the curve-fit.
%               corrP   Correlation matrix of the parameters.
%               Rsq     R-squared cofficient of multiple determination.
%               cvgHst  Convergence history.
%               exitFlag The reason, why the code stopped:
%                           1       convergence in r.h.s. ("JtWdy"),
%                           2       convergence in parameters,
%                           3       convergence in reduced Chi-square,
%                           4       maximum Number of iterations reached
%                                   without convergence
%               message String, one of the above messages.
%               nIter   The number of iterations executed during the fit.
%               nFunEvals The number of function evaluations executed
%                       during the fit.
%
% See also NDBASE.PSO.
%

% Reference:
% Coelho, Alan Anthony. "Optimum Levenberg–Marquardt constant determination
% for nonlinear least-squares." 
% Journal of Applied Crystallography 51.2 (2018): 428-435.


nparams = numel(p0);

inpForm.fname  = {'diff_step' 'lb'            'ub'           'MaxIter' };
inpForm.defval = {1e-7        -inf(1,nparams) inf(1,nparams) 100*nparams};
inpForm.size   = {[1 -1]      [1 nparams]     [1 nparams]    [1 1]};

inpForm.fname  = [inpForm.fname  {'gtol' 'ftol' 'xtol' 'lambda0'}];
inpForm.defval = [inpForm.defval {1e-8    1e-8   1e-8  1e-2}];
inpForm.size   = [inpForm.size   {[1 1]   [1 1]  [1 1] [1 1]}];

inpForm.fname  = [inpForm.fname  {'nu_up', 'nu_dn'}];
inpForm.defval = [inpForm.defval {10       0.3}];
inpForm.size   = [inpForm.size   {[1 1]    [1 1]}];

param = sw_readparam(inpForm, varargin{:});

cost_func_wrap = ndbase.cost_function_wrapper(func, p0, "data", dat, 'lb', param.lb, 'ub', param.ub);
% transform starting values into their unconstrained surrogates.
p0_free = cost_func_wrap.get_free_parameters(p0);

if isempty(dat)
    % minimising scalar - empty resid
    eval_cost_func = @eval_cost_scalar;
    calc_hessian_and_jacobian = @calc_hessian_and_jacobian_scalar;
    ndof = 1;
    diff_step = param.diff_step; % always interpreted as fractional step even if 1 parameter
else
    % minimising sum square residuals
    eval_cost_func = @eval_cost_resids;
    calc_hessian_and_jacobian = @calc_hessian_and_jacobian_resid;
    ndof = numel(dat.x) - cost_func_wrap.get_num_free_parameters() + 1;
    % get absolute diff_step for each parameter
    diff_step = abs(p0_free).*param.diff_step;
    min_step = sqrt(eps);
    diff_step(abs(diff_step) < min_step) = min_step;
end

if isempty(p0_free)
    % All parameters fixed, evaluate cost at initial guess
    % don't use p0 as could contain fixed params outside bounds
    pOpt = cost_func_wrap.get_bound_parameters(p0_free);
    fVal = cost_func_wrap.eval_cost_function(p0_free)/ndof;
    message = 'Parameters are fixed, no optimisation';
    perr = zeros(size(pOpt));
    niter = 0;
    exit_flag = 0;
else
    p = p0_free(:);
    lambda = param.lambda0;
    [cost_val, resids] = eval_cost_func(cost_func_wrap, p);
    [hess, jac] = calc_hessian_and_jacobian_scalar(cost_func_wrap, p, diff_step, cost_val, resids);
    exit_flag = 0;
    for niter = 1:param.MaxIter
        dp = calc_parameter_step(hess, jac, lambda);
        if norm(dp) < param.xtol*(param.xtol + norm(p))
            message = "step size below tolerance xtol";
            exit_flag = 1;
            break
        end
        new_p = p - dp;
        new_cost_val = cost_func_wrap.eval_cost_function(new_p);
        dcost = new_cost_val - cost_val;
        if dcost < 0
            lambda = param.nu_dn*lambda; % decrease step
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
            % limit number times do this?
            lambda = param.nu_up*lambda;
        end
    end
    % collect output
    pOpt = cost_func_wrap.get_bound_parameters(p);
    fVal = cost_val / ndof;
    if exit_flag > 0
        % converged on solution - calculate errors
        cov = pinv(hess) * 2.0 * fVal;
        perr = sqrt(diag(cov));
    else
        message = "Failed to converge in MaxIter";
        perr = zeros(size(pOpt));
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

