function hessian = estimate_hessian(fcost, params, step_size)
% ### Syntax
% 
% `hessian = ndbase.estimate_hessian(func_cost, params)`
%
% `hessian = ndbase.estimate_hessian(func_cost, params, step_size)`
% 
% ### Description
% 
% Function to estimate hessian (curvature) matrix at the best-fit
% parameters given a callable cost function, using finite-difference method.
% The results have been checked against MATLAB's lsqnonlin.
%
% The covariance matrix can be calculated from the hessian
% 
% `cov = inv(hess) * 2 * chisq_red`
%
% where `chisq_red` is the reduced chi-squared at the best fit parameters.
%
% The uncertainty on the best fit parameters is then given by
% `param_errors = sqrt(diag(cov))`
%
% Note this method is more memory efficient than the error estimation in 
% e.g. ndbase.lm3 which uses the Gauss-Newton approx. to the hessian from
% the partial derivative of the fit function with respect to each parameter
% (jacobian of size [ndata x nparams]). However, evaluating the hessian
% explicitly requires more function calls.
% 
% ### Input Arguments
% 
% `fcost`
% : Callable cost function of form `fcost = @(params) ...` taking vector 
%   of fit parameters, `params` and returning a scalar
%
% `params`
% : Vector of best-fit parameters
%
% `step_size` (optional)
% : Scalar double or vector of doubles (same size as params) governing the
%   step used in the finite difference derivatives:
%
%   * If a scalar step_size is provided then this will be interpreted as a
%     fractional step_size i.e. the step used will be step_size*parameter.
%   * If a vector step_size is provided this will denote the absolute step
%     size for each parameter.
%   * If step_size is not provided the step size will be optimised
%     starting from a minimum fractional step_size of sqrt(eps) ~ 1e-8.
%     Note this may require additional function evauations - which may be
%     expensive.
%
% ### Output Arguments
% 
% `hessian` 
% : Hessian matrix from which covariance matrix can be calculated.
%
% ### Examples
%```
% >> % make data
% >> rng default;
% >> fit_fun = @(x, params) params(1)*x + params(2)./x + params(3);
% >> pars = [5,4,3];
% >> x = linspace(0.1,2);
% >> y = fit_fun(x, pars);
% >> e = sqrt(y);
% >> y = y + e.*randn(size(y));
% >>
% >> % fit data
% >> fcost = @(params) sum(((y-fit_fun(x, params))./e).^2);  % chisq
% >> [pars_fit, chisq, ~] = ndbase.simplex([],fcost,pars);
% >>
% >> % estimate errors on fit parameters
% >> hess = ndbase.estimate_hessian(fcost, pars_fit);
% >> cov = inv(hess) * 2.0 * chisq/ (numel(y) - numel(pars_fit));
% >> perr = sqrt(diag(cov)); % errors on best fit parameters
%'''
%
    npar = numel(params);
    optimise_step = false;
    min_step = sqrt(eps);
    max_niter = 16; % default number of iterations in optimising step size
    if nargin > 2
        if numel(step_size)==1
            % interpret as fractional step size
            step_size = params.*step_size;
        elseif numel(step_size) ~= npar
            error("ndbase:estimate_hessian", ...
                  "If step_size is provided it must be a scalar or vector " + ...
                  "with same number of elements as params");
        end
    else
        step_size = params.*min_step;
        optimise_step=true;
    end
    step_size(step_size < min_step) = min_step;
    
    % calculate jacobian at param
    hessian = zeros(npar);
    initial_cost = fcost(params);
    jac_one_step = zeros(1, npar);
    cost_one_step = zeros(1, npar);
    for ipar = 1:npar
        success = false;
        for niter = 1:max_niter
            params(ipar) = params(ipar) + step_size(ipar);
            cost_one_step(ipar) = fcost(params);
            delta_cost = (cost_one_step(ipar) - initial_cost);
            params(ipar) = params(ipar) - step_size(ipar);
            if delta_cost > min_step || ~optimise_step
                success = true;
                break
            else
                step_size(ipar) = 2*step_size(ipar);
            end
        end
        if optimise_step && ~success
            warning("ndbase:estimate_hessian", "Failed to optimise step size for parameter index %d", ipar);
        end
        jac_one_step(ipar) = delta_cost/step_size(ipar);
    end
    % compute hessian
    for irow = 1:npar
        params(irow) = params(irow) + step_size(irow);
        for icol = irow:npar
            params(icol) = params(icol) + step_size(icol);
            cost_two_step = fcost(params);
            jac_two_step = (cost_two_step - cost_one_step(irow))/step_size(icol);
            hessian(irow, icol) = (jac_two_step - jac_one_step(icol))/step_size(irow);
            hessian(icol,irow) = hessian(irow, icol);
            params(icol) = params(icol) - step_size(icol);
        end
        params(irow) = params(irow) - step_size(irow);
    end

end
