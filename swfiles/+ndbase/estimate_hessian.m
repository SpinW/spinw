function [hessian, varargout] = estimate_hessian(fcost, params, options)
% ### Syntax
% 
% `hessian = ndbase.estimate_hessian(func_cost, params)`
%
% `[hessian, cost_val] = ndbase.estimate_hessian(func_cost, params)`
%
% `hessian = ndbase.estimate_hessian(func_cost, params, 'step', step_size)`
%
% `hessian = ndbase.estimate_hessian(func_cost, params, 'ivary', ivary)`
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
% ### Name-Value Pair Arguments
%
% `step` (optional)
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
% `ivary` (optional)
% : indices of parameters which were varied in the fit - the final hessian
%   will be a sqare matrix of size [numel(ivary) x numel(ivary)].
%   If not provided all parameters will be varied.
%
% `niter` (optional)
% : Number of iterations in optimising step size (step doubled in each
%   iteration from minimum value of sqrt(eps))
%
% `cost_tol` (optional)
% : Minimum difference in cost function value for a step size to be
%   considered appropriate. Step size optimisation will be terminated once
%   this condition is met. Default is sqrt(epse) ~ 1e-8.
%
% ### Output Arguments
% 
% `hessian` 
% : Hessian matrix from which covariance matrix can be calculated.
%
% `cost_val` 
% : Value of cost function evaluated at input parameters in params vector.
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
    arguments
        fcost function_handle
        params double
        options.step double = 0
        options.ivary double = 0
        options.niter (1,1) double = 16
        options.cost_tol (1,1) double = sqrt(eps)
    end
    optimise_step = true;
    min_step = sqrt(eps);
    step_size = params.*min_step;
    npar = numel(params);
    ivary = 1:npar;
    if options.step
        optimise_step = false;
        step_size = options.step;
        if numel(step_size)==1
            % interpret as fractional step size
            step_size = params .* step_size;
        end
    end
    if options.ivary
        ivary = options.ivary;
        npar = numel(ivary);
    end
    max_niter = options.niter;
    cost_tol = options.cost_tol;


    if numel(step_size) ~= numel(params)
        error("ndbase:estimate_hessian", ...
              "If step_size is provided it must be a scalar or vector " + ...
              "with same number of elements as params");
    end
    % enforce minimum step size
    step_size(abs(step_size) < min_step) = min_step;
    
    % calculate jacobian at param
    hessian = zeros(npar);
    initial_cost = fcost(params);
    jac_one_step = zeros(1, npar);
    cost_one_step = zeros(1, npar);
    for irow = 1:npar
        ipar = ivary(irow);
        success = false;
        for niter = 1:max_niter
            params(ipar) = params(ipar) + step_size(ipar);
            cost_one_step(irow) = fcost(params);
            delta_cost = (cost_one_step(irow) - initial_cost);
            params(ipar) = params(ipar) - step_size(ipar);
            if abs(delta_cost) > cost_tol || ~optimise_step
                success = true;
                break
            else
                step_size(ipar) = 2*step_size(ipar);
            end
        end
        if optimise_step && ~success
            warning("ndbase:estimate_hessian", "Failed to optimise step size for parameter index %d", ipar);
        end
        jac_one_step(irow) = delta_cost/step_size(ipar);
    end
    % compute hessian
    for irow = 1:npar
        ipar_row = ivary(irow);
        params(ipar_row) = params(ipar_row) + step_size(ipar_row);
        for icol = irow:npar
            ipar_col = ivary(icol);
            params(ipar_col) = params(ipar_col) + step_size(ipar_col);
            cost_two_step = fcost(params);
            jac_two_step = (cost_two_step - cost_one_step(irow))/step_size(ipar_col);
            hessian(irow, icol) = (jac_two_step - jac_one_step(icol))/step_size(ipar_row);
            hessian(icol,irow) = hessian(irow, icol);
            params(ipar_col) = params(ipar_col) - step_size(ipar_col);
        end
        params(ipar_row) = params(ipar_row) - step_size(ipar_row);
    end
    varargout = {initial_cost};
end
