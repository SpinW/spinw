classdef sw_fitpowder < handle & matlab.mixin.SetGet
    % class to fit powders
    %
    % [spinw.copy], [spinw.struct], [Comparing handle and value classes](https://www.mathworks.com/help/matlab/matlab_oop/comparing-handle-and-value-classes.html)
    %

    properties (SetObservable)
       swobj;
       y
       e
       ebin_cens
       modQ_cens
       nQ = 10  % number |Q| bins to calc. in integration limits of 1D cut
       ndim
       ncuts = 0;
       sw_instrument_args = struct()
       powspec_args = struct('nRand', 1e3, 'T', 0, 'formfact', true, 'hermit', false,  'imagChk', false, 'fibo', true, 'binType', 'cbin', 'component', 'Sperp')
       fit_func
       cost_function = "Rsq";  % "Rsq" or "chisq"
       background_strategy = "planar" % "planar" or "independent" (1D only - fbg = @(en, p1, p2, ..., pN)
       fbg = @(en, modQ, slope_en, slope_modQ, intercept) slope_en*en(:) + slope_modQ*modQ(:)' + intercept
       % parameters
       nparams_model
       nparams_bg
       params
       bounds
    end
        
    methods       
        function obj = sw_fitpowder(swobj, data, fit_func, model_params, varargin)
            % constructor
            obj.swobj = swobj;
            obj.fit_func = fit_func;
            bg_params = [];
            scale = 1;
            % set background startegy and func
            for ikey = 1:2:numel(varargin)
                switch varargin{ikey}
                    case 'backgroundStrategy'
                        obj.background_strategy = varargin{ikey+1};
                    case 'initialBackgroundParameters'
                        bg_params = varargin{ikey+1};
                    case 'scale'
                        scale = varargin{ikey+1};
                end
            end
            % initialise parameters
            obj.add_data(data)
            obj.initialise_parameters_and_bounds(model_params, bg_params, scale)
        end

        function initialise_parameters_and_bounds(obj, model_params, bg_params, scale)
            obj.nparams_model = numel(model_params);
            obj.nparams_bg = obj.get_nparams_in_background_func();
            if isempty(bg_params)
                % zero intialise
                if obj.background_strategy == "independent"
                    nparams_bg_total = obj.ncuts * obj.nparams_bg;
                else
                    nparams_bg_total = obj.nparams_bg;
                end
                bg_params = zeros(nparams_bg_total, 1);
            end
            obj.params = [model_params(:); bg_params(:); scale];
            obj.bounds = [-inf, inf].*ones(numel(obj.params),1);
            obj.bounds(end,1) = 0; % set lower bound of scale to be 0
        end

        function fix_parameter(obj, iparams, varargin)
            for ii = 1:numel(iparams)
                iparam = obj.get_parameter_index(iparams(ii), varargin{:});
                obj.bounds(iparam, :) = obj.params(iparam);
            end
        end

        function set_parameter(obj, iparams, values, varargin)
            for ii = 1:numel(iparams)
                iparam = obj.get_parameter_index(iparams(ii), varargin{:});
                obj.params(iparam) = values(ii);
            end
        end

        function set_parameter_bounds(obj, iparams, varargin)
            lb = [];
            ub = [];
            for ikey = 1:2:numel(varargin)
                switch varargin{ikey}
                    case 'lb'
                        lb = varargin{ikey+1};
                    case 'ub'
                        ub = varargin{ikey+1};
                end
            end
            for ii = 1:numel(iparams)
                iparam = obj.get_parameter_index(iparams(ii), varargin{:});
                if ~isempty(lb)
                    obj.bounds(iparam, 1) = lb(ii);
                end
                if ~isempty(ub)
                    obj.bounds(iparam, 2) = ub(ii);
                end
            end
        end

        function add_data(obj, data)
            if numel(data) == 1 && numel(data.x) == 2
                % 2D
                obj.ndim = 2;
                obj.y = data.y';  % nE x n|Q|
                obj.e = data.e';
                obj.ebin_cens = data.x{1};
                obj.modQ_cens = data.x{2};
            else
                % 1D
                obj.ndim = 1;
                obj.ebin_cens = data(1).x;
                obj.ncuts = numel(data);
                for cut = data
                    obj.y = [obj.y cut.y(:)];
                    obj.e = [obj.e cut.e(:)];
                    obj.modQ_cens = [obj.modQ_cens, linspace(cut.qmin, cut.qmax, obj.nQ)];  % does this play nicely with sw_instrument Q convolution?
                end
            end
            obj.powspec_args.Evect = obj.ebin_cens(:)';  % row vect
        end

        function crop_energy_range(obj, emin, emax)
            % crop data
            ikeep = obj.ebin_cens >= emin & obj.ebin_cens <= emax;
            obj.ebin_cens = obj.ebin_cens(ikeep);
            obj.y = obj.y(ikeep, :);
            obj.e = obj.e(ikeep, :);
        end

        function bg = calc_background(obj, params)
            % add background
            istart = obj.nparams_model + 1;
            if obj.background_strategy == "planar"
                bg_params = num2cell(params(istart:istart+obj.nparams_bg-1));
                bg = obj.fbg(obj.ebin_cens, obj.modQ_cens, bg_params{:});
            elseif obj.ndim == 1
                % add energy dependent background to each cut
                bg = zeros(size(obj.y));
                for icut = 1:obj.ncuts
                    iend = istart + +obj.nparams_bg - 1;
                    bg_params = num2cell(params(istart:iend));
                    bg(:,icut) = obj.fbg(obj.ebin_cens(:), bg_params{:});
                    istart = iend + 1;
                end
            end
        end

        function [ycalc, bg] = calc_spinwave_spec(obj, params)
            % set spinw model interactions
            obj.fit_func(obj.swobj, params(1:obj.nparams_model)');
            % simulate powder spectrum
            spec = obj.swobj.powspec(obj.modQ_cens, obj.powspec_args);
            spec = sw_instrument(spec, obj.sw_instrument_args);
            % scale
            ycalc = params(end)*spec.swConv;
            % add background
            bg = calc_background(obj, params);
            ycalc = ycalc + bg;
            if obj.ndim == 1
                % integrate nQ |Q| points for each cut
                ycalc = obj.rebin_powspec_to_1D_cuts(ycalc);
            end
        end

        function resid_sq = calc_resid_sq(obj, params)
            % evaluate fit function
            [ycalc, ~] = calc_spinwave_spec(obj, params);
            resid_sq = (obj.y - ycalc).^2;
            if obj.cost_function == "chisq"
                resid_sq = resid_sq./(obj.e.^2);
            end
            % exclude nans in both ycalc and input data
            resid_sq = resid_sq(~isnan(resid_sq));
        end

        function result = fit(obj, options)           
            % call fit
            [p, cost, resid_sq, exitflag, ~,~,jacobian] = lsqnonlin(@obj.calc_resid_sq,obj.params,obj.bounds(:,1),obj.bounds(:,2), options);
            result.p = p;
            result.exitflag = exitflag;
            % calculate errors
            ndof = numel(resid_sq) - numel(p);
            result.reduced_cost = (cost/ndof);
            cov = result.reduced_cost * inv(jacobian' * jacobian);
            result.perr = diag(cov);            
        end

        function estimate_scale_factor(obj)
            % set scale factor to 1
            params = obj.params;
            params(end) = 1;
            [ycalc, bg] = calc_spinwave_spec(obj, params);
            if obj.ndim == 1
                % integrate nQ |Q| points for each cut
                bg = obj.rebin_powspec_to_1D_cuts(bg);
            end
            scale = max(obj.y - bg, [], "all")/max(ycalc, [], "all");
            obj.params(end) = scale;
        end

    end

    % private
    methods (Hidden=true, Access = private)
        function ycalc = rebin_powspec_to_1D_cuts(obj, ycalc)
            % sum up successive nQ points along |Q| axis (dim=2)
            ycalc = reshape(ycalc, size(ycalc,1), obj.nQ, []);
            ycalc = squeeze(sum(ycalc, 2, 'omitnan'));
        end

        function nbg_pars = get_nparams_in_background_func(obj)
            % get background parameters
            if obj.background_strategy == "planar"
                nbg_pars = nargin(obj.fbg) - 2; % dependent variables energy, modQ
            else
                nbg_pars = nargin(obj.fbg) - 1; % dependent variable energy
            end
        end

        function iparam = get_parameter_index(obj, iparam, varargin)
            icut = 0;
            is_background_param = false;
            for ikey = 1:2:numel(varargin)
                switch varargin{ikey}
                    case 'background'
                        is_background_param = varargin{ikey+1};
                    case 'cutIndex'
                        if obj.ndim == 1 && obj.background_strategy=="independent"
                            icut = varargin{ikey+1};
                            if icut < 1 || icut > obj.ncuts
                                error("Invalid CutIndex supplied.")
                            end
                        else
                            error("Can only supply cutIndex when fitting 1D cuts with independent backgrounds")
                        end
                end
            end
            % find index in vector of all parameters
            if is_background_param
                iparam = iparam + obj.nparams_model;
                if obj.background_strategy == "independent"
                    nbg_pars = obj.get_nparams_in_background_func();
                    iparam = iparam + (icut-1)*nbg_pars;
                end
            end
        end
    end
end