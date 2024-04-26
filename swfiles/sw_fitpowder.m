classdef sw_fitpowder < handle & matlab.mixin.SetGet
    % class to fit powders
    %
    % [spinw.copy], [spinw.struct], [Comparing handle and value classes](https://www.mathworks.com/help/matlab/matlab_oop/comparing-handle-and-value-classes.html)
    %

    properties (SetObservable)
        % data
       swobj;
       y
       e
       ebin_cens
       modQ_cens
       nQ = 10  % number |Q| bins to calc. in integration limits of 1D cut
       ndim
       ncuts = 0;
       % spinw funtion inputs
       sw_instrument_args = struct()
       powspec_args = struct('nRand', 1e3, 'T', 0, 'formfact', true, 'hermit', false,  'imagChk', false, 'fibo', true, 'binType', 'cbin', 'component', 'Sperp')
       % functions
       fit_func
       cost_function = "Rsq";  % "Rsq" or "chisq"
       background_strategy = "planar" % "planar" or "independent" (1D only - fbg = @(en, p1, p2, ..., pN)
       fbg
       fbg_planar = @(en, modQ, slope_en, slope_modQ, intercept) slope_en*en(:) + slope_modQ*modQ(:)' + intercept;
       fbg_indep = @(en, slope_en, intercept) slope_en*en(:) + intercept;
       % fit and parameters
       optimizer = @ndbase.simplex
       nparams_model
       nparams_bg
       params
       bounds
       % cache calculation
       do_cache = true
       ycalc_cached = []
       model_params_cached = []
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
            obj.add_data(data)
            if obj.background_strategy == "planar"
                obj.fbg = obj.fbg_planar;
            else
                obj.fbg = obj.fbg_indep;
            end
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

        function fix_model_parameters(obj, iparams)
            if any(iparams > obj.nparams_model)
                error('sw_fitpowder', 'Parameter indices supplied must be within number of model parameters');
            end
            for iparam = iparams
               obj.bounds(iparam, :) = obj.params(iparam);
            end
        end

        function fix_bg_parameters(obj, iparams_bg, icuts)
            if any(iparams_bg > obj.nparams_bg)
                error('sw_fitpowder', 'Parameter indices supplied must be within number of background parameters');
            end
            if nargin < 3
                icuts = 0;
            end
            for iparam = obj.get_index_of_background_parameters(iparams_bg, icuts)
               obj.bounds(iparam, :) = obj.params(iparam);
            end
        end

        function fix_scale(obj)
            obj.bounds(end, :) = obj.params(end);
        end

        function set_model_parameters(obj, iparams, values)
            if any(iparams > obj.nparams_model)
                error('sw_fitpowder', 'Parameter indices supplied must be within number of model parameters');
            end
            for ival = 1:numel(values)
               obj.params(iparams(ival)) = values(ival);
            end
        end

        function set_bg_parameters(obj, iparams_bg, values, icuts)
            if any(iparams_bg > obj.nparams_bg)
                error('sw_fitpowder', 'Parameter indices supplied must be within number of background parameters');
            end
            if nargin < 4
                icuts = 0;
            end
            for ival = 1:numel(values)
                iparams = obj.get_index_of_background_parameters(iparams_bg(ival), icuts);
                obj.params(iparams) = values(ival);
            end
        end

        function set_scale(obj, scale)
            obj.params(end) = scale;
        end

        function set_model_parameter_bounds(obj, iparams, lb, ub)
            if any(iparams > obj.nparams_model)
                error('sw_fitpowder', 'Parameter indices supplied must be within number of model parameters');
            end
            for ibnd = 1:numel(iparams)
                if ~isempty(lb)
                    obj.bounds(iparams(ibnd), 1) = lb(ibnd);
                end
                if ~isempty(ub)
                    obj.bounds(iparams(ibnd), 2) = ub(ibnd);
                end
            end
        end

        function set_bg_parameter_bounds(obj, iparams_bg, lb, ub, icuts)
            if any(iparams_bg > obj.nparams_bg)
                error('sw_fitpowder', 'Parameter indices supplied must be within number of background parameters');
            end
            if nargin < 5
                icuts = 0;
            end
            for ibnd = 1:numel(iparams_bg)
                iparams = obj.get_index_of_background_parameters(iparams_bg(ibnd), icuts);
                if ~isempty(lb)
                    obj.bounds(iparams, 1) = lb(ibnd);
                end
                if ~isempty(ub)
                    obj.bounds(iparams, 2) = ub(ibnd);
                end
            end
        end

        function set_scale_bounds(obj, lb, ub)
            if ~isempty(lb)
                obj.bounds(end, 1) = lb;
            end
            if ~isempty(ub)
                obj.bounds(end, 2) = ub;
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
            obj.powspec_args.Evect = obj.ebin_cens(:)';
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

        function set_caching(obj, do_cache)
            obj.do_cache = do_cache;
            if ~obj.do_cache
                obj.clear_cache()
            end
        end

        function clear_cache(obj)
            obj.ycalc_cached = [];
            obj.model_params_cached = [];
        end

        function [ycalc, bg] = calc_spinwave_spec(obj, params)
            model_params = reshape(params(1:obj.nparams_model), 1, []);
            if obj.do_cache && ~isempty(obj.model_params_cached) && all(abs(model_params - obj.model_params_cached) < 1e-10)
                ycalc = obj.ycalc_cached;  % do not repeat calc
            else
                % set spinw model interactions
                % pass params as row-vector as expected by matparser
                obj.fit_func(obj.swobj, model_params);
                % simulate powder spectrum
                spec = obj.swobj.powspec(obj.modQ_cens, obj.powspec_args);
                spec = sw_instrument(spec, obj.sw_instrument_args);
                ycalc = spec.swConv;
                % cache if required
                if obj.do_cache
                    obj.ycalc_cached = spec.swConv;
                    obj.model_params_cached = model_params;
                end
            end
            % scale
            ycalc = params(end)*ycalc;
            % calc background
            bg = calc_background(obj, params);
            if obj.background_strategy == "planar"
                ycalc = ycalc + bg;  % add planar bg before any rebinning
                if obj.ndim == 1
                    % integrate nQ |Q| points for each cut
                    ycalc = obj.rebin_powspec_to_1D_cuts(ycalc);
                end
            else
                % integrate nQ |Q| points for each cut
                ycalc = obj.rebin_powspec_to_1D_cuts(ycalc);
                % add background to individual cuts after rebin
                ycalc = ycalc + bg;
            end

        end

        function resid_sq_sum = calc_cost_func(obj, params)
            % evaluate fit function
            [ycalc, ~] = calc_spinwave_spec(obj, params);
            resid = (obj.y - ycalc);
            if obj.cost_function == "chisq"
                resid = resid./(obj.e);
            end
            % exclude nans in both ycalc and input data
            ikeep = ~isnan(resid);
            resid_sq_sum = resid(ikeep)'*resid(ikeep);
        end

        function result = fit(obj, varargin) 
            % setup cell for output of ndbase optimizer/minimizer
            result = cell(1,nargout(obj.optimizer));
            % pass params and bounds as rows for ndbase optimisers
            [result{:}] = obj.optimizer([], @obj.calc_cost_func, obj.params(:)', ...
                                       'lb', obj.bounds(:,1)', ...
                                       'ub', obj.bounds(:,2)', ...
                                       varargin{:});
        end

        function estimate_scale_factor(obj)
            % set scale factor to 1
            params = obj.params;
            params(end) = 1;
            [ycalc, bg] = calc_spinwave_spec(obj, params);
            scale = max(obj.y - mean(bg(:)), [], "all")/max(ycalc, [], "all");
            obj.params(end) = scale;
        end

        function plot_result(obj, params, varargin)
            [ycalc, ~] = obj.calc_spinwave_spec(params);
            if obj.ndim == 1
                obj.plot_1d_cuts_on_data(ycalc, varargin{:})
            else
                obj.plot_2d_contour_on_data(ycalc, varargin{:})
            end
        end

        function plot_1d_cuts_on_data(obj, ycalc, varargin)
            figure("color","white");
            modQs = mean(reshape(obj.modQ_cens, [], obj.ncuts), 1);
            for icut = 1:obj.ncuts
                ax =  subplot(1, obj.ncuts, icut);
                hold on; box on;
                plot(ax, obj.ebin_cens, obj.y(:,icut), 'ok');
                plot(ax, obj.ebin_cens, ycalc(:,icut), '-r', varargin{:});
                xlim(ax, [0, 12]);
                ylim(ax, [-2, 30]);
                xlabel(ax, 'Energy (meV)')
                ylabel(ax, 'Intensity');
                title(ax, num2str(modQs(icut), 2) + " $\AA^{-1}$", 'interpreter','latex')
            end
        end

        function plot_2d_contour_on_data(obj, ycalc, varargin)
            % varargin passed to contour
            figure("color","white");
            ax = subplot(1,1,1);
            box on; hold on;
            h = imagesc(ax, obj.modQ_cens, obj.ebin_cens, obj.y);
            h.AlphaData = obj.y > 0;  % make empty bins transparent
            contour(ax, obj.modQ_cens, obj.ebin_cens, ycalc, varargin{:});
            cbar = colorbar(ax);
            cbar.Label.String = "Intensity";
            xlabel(ax, "$\left|Q\right| (\AA^{-1})$", 'interpreter','latex');
            ylabel(ax, "Energy (meV)");
            legend('Calculated');
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

        function iparams = get_index_of_background_parameters(obj, iparams_bg, icuts)
            if any(icuts == 0) && obj.background_strategy == "independent"
                icuts = 1:obj.ncuts;  % apply to all cuts
            end
            % find index in vector of all parameters
            iparams = [];
            for icut = icuts
                iparams = [iparams,  iparams_bg + obj.nparams_model + (icut-1)*obj.nparams_bg];
            end
        end
    end
end