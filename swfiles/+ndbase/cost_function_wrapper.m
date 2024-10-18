classdef cost_function_wrapper < handle & matlab.mixin.SetGet
% ### Syntax
% 
% `param = fit_parameter(value, lb, ub)`
% 
% ### Description
% 
% Class for evaluating cost function given data and parameters.
% Optionally the parameters can be bound in which case the class will
% perform a transformation to convert the constrained optimization problem
% into an un-constrained problem, using the formulation devised
% (and documented) for MINUIT [1] and also used in lmfit [2].
%
% [1] https://root.cern/root/htmldoc/guides/minuit2/Minuit2.pdf#section.2.3
% [2] https://lmfit.github.io/lmfit-py/bounds.html
% 
% ### Input Arguments
% 
% `func`
% : Function handle or char with one of the following definition:
%   * `R2 = func(p)`        if `dat` is empty,
%   * `y  = func(x,p)`      if `dat` is a struct.
%   Here `x` is a vector of $N$ independent variables, `p` are the
%   $M$ parameters to be optimized and `y` is the simulated model, `R2`
%   is the value to minimize.
%
% `parameters`
% :  Vector of doubles
%
% ### Name-Value Pair Arguments
%
% `data`
% : Either empty or contains data to be fitted stored in a structure with
%   fields:
%   * `dat.x`   vector of $N$ independent variables,
%   * `dat.y`   vector of $N$ data values to be fitted,
%   * `dat.e`   vector of $N$ standard deviation (positive numbers)
%               used to weight the fit. If zero or missing
%               an unweighted fit will be performed.
%
% `lb`
% : Optional vector of doubles corresponding to the lower bound of the 
%   parameters. Empty vector [] or vector of -inf interpreted as no lower 
%   bound.
% 
% `ub`
% : Optional vector of doubles corresponding to the upper bound of the 
%   parameters. Empty vector [] or vector of inf interpreted as no upper 
%   bound.
%
% ### Examples
    properties (SetObservable)
        % data
        cost_func
        free_to_bound_funcs
        bound_to_free_funcs
        ifixed
        ifree
        pars_fixed
    end

    methods
        function obj = cost_function_wrapper(fhandle, params, options)
            arguments
                fhandle {isFunctionHandleOrChar}
                params double
                options.lb double = []
                options.ub double = []
                options.data struct = struct()
            end
            if ischar(fhandle)
                fhandle = str2func(fhandle); % convert to fuction handle
            end
            if isempty(fieldnames(options.data))
                % fhandle calculates cost_val
                obj.cost_func = fhandle;
            else
                % fhandle calculates fit/curve function
                calc_resdid = @(p) fhandle(options.data.x(:), p) - options.data.y(:);
                if ~isfield(options.data,'e') || isempty(options.data.e) || ~any(options.data.e(:))
                    warning("ndbase:cost_function_wrapper:InvalidWeights",...
                        "Invalid weights provided - unweighted residuals will be used.")
                    obj.cost_func = @(p) sum(calc_resdid(p).^2);
                else
                    obj.cost_func = @(p) sum((calc_resdid(p)./options.data.e(:)).^2);
                end
            end

            % validate size of bounds
            lb = options.lb;
            ub = options.ub;
            if ~isempty(lb) && numel(lb) ~= numel(params)
                error("ndbase:cost_function_wrapper:WrongInput", ...
                  "Lower bounds must be empty or have same size as parameter vector.");
            end
            if ~isempty(ub) && numel(ub) ~= numel(params)
                error("ndbase:cost_function_wrapper:WrongInput", ...
                  "Upper bounds must be empty or have same size as parameter vector.");
            end
            if ~isempty(lb) && ~isempty(ub) && any(ub<lb)
                error("ndbase:cost_function_wrapper:WrongInput", ...
                    "Upper bounds have to be larger than the lower bounds.");
            end
            % init bound parameters
            obj.init_bound_parameter_transforms(params, lb, ub)
        end

        function init_bound_parameter_transforms(obj, pars, lb, ub)
            % Note free parameters to be used externally in the
            % optimisation (note in lmfit [2] pfree is called p_internal).
            % Bound parameters are the original parameters
            % pass into the constructor (that may or may not be bound or
            % fixed).
            obj.free_to_bound_funcs = cell(size(pars));
            obj.bound_to_free_funcs = cell(size(pars));
            obj.ifixed = [];
            ipars = 1:numel(pars); % used later
            for ipar = ipars
                has_lb = ~isempty(lb) && lb(ipar) > -inf;
                has_ub = ~isempty(ub) && ub(ipar) < inf;
                if has_lb && has_ub
                    % both bounds specified and parameter not fixed
                    obj.free_to_bound_funcs{ipar} = @(p) obj.free_to_bound_has_lb_and_ub(p, lb(ipar), ub(ipar));
                    obj.bound_to_free_funcs{ipar} = @(p) obj.bound_to_free_has_lb_and_ub(p, lb(ipar), ub(ipar));
                    % check if fixed
                    if abs(ub(ipar) - lb(ipar)) < max(abs(ub(ipar)), 1)*1e-10
                        obj.ifixed = [obj.ifixed, ipar];
                        if  pars(ipar) < lb(ipar)
                             pars(ipar) = lb(ipar);
                        elseif  pars(ipar) > ub(ipar)
                             pars(ipar) = ub(ipar);
                        end
                        obj.pars_fixed = [obj.pars_fixed, pars(ipar)];
                    end
                elseif has_lb
                    obj.free_to_bound_funcs{ipar} = @(p) obj.free_to_bound_has_lb(p, lb(ipar));
                    obj.bound_to_free_funcs{ipar} = @(p) obj.bound_to_free_has_lb(p, lb(ipar));
                elseif has_ub
                    obj.free_to_bound_funcs{ipar} = @(p) obj.free_to_bound_has_ub(p, ub(ipar));
                    obj.bound_to_free_funcs{ipar} = @(p) obj.bound_to_free_has_ub(p, ub(ipar));
                else
                    obj.free_to_bound_funcs{ipar} = @(p) p;
                    obj.bound_to_free_funcs{ipar} = @(p) p;
                end
            end
            % get index of free parameters
            obj.ifree = ipars(~ismember(1:numel(pars), obj.ifixed));
        end

        function pars_bound = get_bound_parameters(obj, pars)
            % pars vector will not include fixed parameters
            pars_bound = zeros(size(obj.free_to_bound_funcs));
            for ipar_free = 1:numel(pars)
                ipar_bound = obj.ifree(ipar_free);
                pars_bound(ipar_bound) = obj.free_to_bound_funcs{ipar_bound}(pars(ipar_free));
            end
            % add in fixed parameter values
            pars_bound(obj.ifixed) = obj.pars_fixed;
        end

        function pars = get_free_parameters(obj, pars_bound)
            pars = zeros(size(pars_bound)); % to preserve par vector shape
            for ipar = obj.ifree
                pars(ipar) = obj.bound_to_free_funcs{ipar}(pars_bound(ipar));
            end
            pars = pars(obj.ifree);
        end

        function cost_val = eval_cost_function(obj, pars_excl_fixed)
            pars = obj.get_bound_parameters(pars_excl_fixed);
            cost_val = obj.cost_func(pars);
        end

        function nfree = get_num_free_parameters(obj)
            nfree = numel(obj.ifree);
        end
    end
    % private
    methods (Static=true, Hidden=true, Access = private)
        function par_bound = free_to_bound_has_lb(par, lb)
            par_bound = lb - 1 + sqrt(par^2 + 1);
        end
        function par = bound_to_free_has_lb(par_bound, lb)
            if par_bound < lb
                par_bound = lb;
            end
            par = sqrt((par_bound - lb + 1)^2 - 1);
        end
        function par_bound = free_to_bound_has_ub(par, ub)
           par_bound = ub + 1 - sqrt(par^2 + 1);
        end
        function par = bound_to_free_has_ub(par_bound, ub)
            if par_bound > ub
                par_bound = ub;
            end
            par = sqrt((ub - par_bound + 1)^2 - 1);
        end
        function par_bound = free_to_bound_has_lb_and_ub(par, lb, ub)
            par_bound = lb + (sin(par) + 1)*(ub-lb)/2;
        end
        function par = bound_to_free_has_lb_and_ub(par_bound, lb, ub)
            if par_bound < lb
                par_bound = lb;
            elseif par_bound > ub
                par_bound = ub;
            end
            par = asin((2*(par_bound - lb)/(ub-lb)) - 1);
        end
    end
end


function isFunctionHandleOrChar(func)
    assert(isa(func, 'function_handle') || ischar(func));
end