classdef sw_fitpowder < handle & matlab.mixin.SetGet
% Class to fit powder averaged spectra to inelastic neutron scattering data
% 
% ### Syntax
% 
% `fitpow = sw_fitpowder(swobj, data, fit_func, model_params, Name,Value)`
% 
% ### Description
% 
% Fits powder averaged spectra to constant-|Q| cuts or 2D |Q| vs en slices
% of inelastic neutron scattering data accounting for the instrument 
% resolution  
% 
% ### Input Arguments
% 
% `swobj`
% : spinwave object with magnetic structure defined
%
% `data`
% : Possible inputs depend on dimensionality of the data:
%   * 2D data 
%     Either a struct containing fields `x`, `y` and `e` or a 2D HORACE 
%     object (sqw or d2d)
%     where `y` is a matrix of intensities with shape (N|Q|, NEnergy)
%           `e` is a matrix of errorsbars with shape (N|Q|, NEnergy)
%           `x` is a cell array of size (1,2): 
%               x{1} contains a vector of energy bin centers
%               x{2} contains a vector of |Q| bin centers.
%   * Vector of 1d datasets at constant |Q|
%     Either a struct containing fields `x`, `y` and `e`  `qmin` and `qmax`
%     or a 1D HORACE object (sqw or d1d)
%     where `y` is a vector of intensities with shape (1, NEnergy)
%           `e` is a vector of errorsbars with shape (1, NEnergy)
%           `x` is a vector of energy bin centers
%               x{2} contains a vector of |Q| bin centers.
%           `qmin` is a scalar denoting the lower Q value of the cut
%           `qmax` is a scalar denoting the upper Q value of the cut
% 
% `fit_func`
% : Function handle changing the interactions in the spinwave model.
%   For example if the spinw model had matrices `J1`, `J2` and `D` then a
%   `fit_func` could be
%   ```
%   fit_func =  @(obj, p) matparser(obj, 'param', p, 'mat', 
%                                   {'J1', 'J2', 'D(3,3)'}, 'init', true);
%   ```
%
% `model_params`
% : Vector of initial paramers to pass to fit_func
%
% `background_strategy` (optional)
% : A string determining the type of background:
%   * `planar` (default) - 2D planar background in |Q| and energy transfer
%   * `independent` - 1D linear background as function of energy transfer
%
% `nQ` (optional)
% : Scalar int correpsonding to number of Q-points to subdivide cuts into 
%   for averaging
%
% ### Output Arguments
% 
% `'result'` 
% : cell array containing output of sw_fitpowder.optimizer
%   For ndbase optimizers in the spinw package then the reuslt can be
%   unpacked as follows
%   ```
%   [fitted_params, cost_val, stat] = result{:}
%   ```
%   See docs for e.g. ndbase.simplex (default) for details
% 
% ### Examples
% 
% ```
% >> % init spinw object
% >> J1 = -0.05;
% >> J2 = 0.3;
% >> D = 0.05;
% >> mnf2 = spinw;
% >> mnf2.genlattice('lat_const', [4.87 4.87 3.31], 'angle', [90 90 90]*pi/180, 'sym', 'P 42/m n m');
% >> mnf2.addatom('r', [0 0 0], 'S', 2.5, 'label', 'MMn2', 'color', 'b')
% >> mnf2.gencoupling('maxDistance', 5)
% >> mnf2.addmatrix('label', 'J1', 'value', J1, 'color', 'red');
% >> mnf2.addmatrix('label', 'J2', 'value', J2, 'color', 'green');
% >> mnf2.addcoupling('mat', 'J1', 'bond', 1)
% >> mnf2.addcoupling('mat', 'J2', 'bond', 2)
% >> mnf2.addmatrix('label', 'D', 'value', diag([0 0 D]), 'color', 'black');
% >> mnf2.addaniso('D')
% >> mnf2.genmagstr('mode', 'direct', 'S', [0 0; 0 0; 1 -1])
% >> 
% >> % define fit_func
% >> fit_func =  @(obj, p) matparser(obj, 'param', p, 'mat', {'J1', 'J2', 'D(3,3)'}, 'init', true);
% >> 
% >> % define resolution (from PyChop)
% >> Ei = 20;
% >> eres = @(en) 512.17*sqrt((Ei-en).^3 .* ( 8.26326e-10*(0.169+0.4*(Ei./(Ei-en)).^1.5).^2 + 2.81618e-11*(1.169+0.4*(Ei./(Ei-en)).^1.5).^2));
% >> % Q-resolution (parameters for MARI)
% >> e2k = @(en) sqrt( en .* (2*1.67492728e-27*1.60217653e-22) )./1.05457168e-34./1e10;
% >> L1 = 11.8;  % Moderator to Sample distance in m
% >> L2 = 2.5;   % Sample to detector distance in m
% >> ws = 0.05;  % Width of sample in m
% >> wm = 0.12;  % Width of moderator in m
% >> wd = 0.025; % Width of detector in m
% >> ki = e2k(Ei);
% >> a1 = ws/L1; % Angular width of sample seen from moderator
% >> a2 = wm/L1; % Angular width of moderator seen from sample
% >> a3 = wd/L2; % Angular width of detector seen from sample
% >> a4 = ws/L2; % Angular width of sample seen from detector
% >> dQ = 2.35 * sqrt( (ki*a1)^2/12 + (ki*a2)^2/12 + (ki*a3)^2/12 + (ki*a4)^2/12 );
% >> 
% >> % fit powder
% >> backgroundStrategy = 'planar';  % 'planar' or 'independent'
% >> nQ = 10;   % Number of Q-points to subdivide cuts into for averaging
% >>
% >> fitpow = sw_fitpowder(mnf2, mnf2dat,  fit_func, [J1 J2 D], backgroundStrategy, nQ);
% >> fitpow.crop_energy_range(2.0, 8.0);
% >> fitpow.powspec_args.dE = eres;
% >> fitpow.sw_instrument_args = struct('dQ', dQ, 'ThetaMin', 3.5, 'Ei', Ei);
% >> fitpow.estimate_constant_background();
% >> fitpow.estimate_scale_factor();
% >> fitpow.set_model_parameter_bounds(1:3, [-1 0 -0.2], [0 1 0.2]) % Force J1 to be ferromagnetic to agree with structure
% >> fitpow.set_bg_parameter_bounds(3, 0.0, []) % Set lb of constant bg = 0
% >> fitpow.fix_bg_parameters(1:2); % fix slopes of background to 0
% >> fitpow.cost_function = "Rsq";  % or "chisq"
% >> fitpow.optimizer = @ndbase.simplex;
% >> result = fitpow.fit('MaxIter', 1);  % passes varargin to optimizer
% >> 
% >> [pfit,cost_val,stat] = result{:};
% >> fitpow.plot_result(pfit, 26, 'EdgeAlpha', 0.9, 'LineWidth', 2)
% ```
%
% [spinw.spinwave] \| [sw_fitspec]
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
       powspec_args = struct('nRand', 1e3, 'T', 0, 'formfact', true, ...
                            'hermit', false,  'imagChk', false, ...
                            'fibo', true, 'binType', 'cbin', ...
                            'component', 'Sperp', 'neutron_output', true, ...
                            'fastmode', true);
       % functions
       fit_func
       cost_function = "Rsq";  % "Rsq" or "chisq"
       background_strategy = "planar" % "planar" or "independent" (1D only - fbg = @(en, p1, p2, ..., pN)
       fbg
       % fit and parameters
       optimizer = @ndbase.simplex
       nparams_model
       nparams_bg
       params
       bounds
       liveplot_interval = 0
    end

    properties (SetAccess = private)
       fbg_planar = @(en, modQ, slope_en, slope_modQ, intercept) slope_en*en(:) + slope_modQ*modQ(:)' + intercept;
       fbg_indep = @(en, slope_en, intercept) slope_en*en(:) + intercept;
       do_cache = true
       ycalc_cached = []
       model_params_cached = []
       liveplot_counter = 0
    end
        
    methods       
        function obj = sw_fitpowder(swobj, data, fit_func, model_params, background_strategy, nQ)
            % constructor
            obj.swobj = swobj;
            obj.fit_func = fit_func;
            if nargin == 6
                obj.nQ = nQ; % set this before add data
            end
            obj.add_data(data)
            if nargin < 5
                obj.background_strategy = "planar";
            else
                obj.background_strategy = background_strategy;
            end
            if obj.background_strategy == "planar"
                obj.fbg = obj.fbg_planar;
            else
                obj.fbg = obj.fbg_indep;
            end
            obj.initialise_parameters_and_bounds(model_params)
        end

        function initialise_parameters_and_bounds(obj, model_params)
            obj.nparams_model = numel(model_params);
            obj.nparams_bg = obj.get_nparams_in_background_func();
            % zero intialise background parameters
            if obj.background_strategy == "independent"
                nparams_bg_total = obj.ncuts * obj.nparams_bg;
            else
                nparams_bg_total = obj.nparams_bg;
            end
            bg_params = zeros(nparams_bg_total, 1);
            % last parameter is scale factor (default to 1)
            obj.params = [model_params(:); bg_params(:); 1];
            obj.bounds = [-inf, inf].*ones(numel(obj.params), 1);
            obj.bounds(end,1) = 0; % set lower bound of scale to be 0
        end

        function fix_model_parameters(obj, iparams)
            if any(iparams > obj.nparams_model)
                error('sw_fitpowder', 'Parameter indices supplied must be within number of model parameters');
            end
            obj.fix_parameters( iparams)
        end

        function fix_bg_parameters(obj, iparams_bg, icuts)
            if any(iparams_bg > obj.nparams_bg)
                error('sw_fitpowder', 'Parameter indices supplied must be within number of background parameters');
            end
            if nargin < 3
                icuts = 0;
            end
            obj.fix_parameters(obj.get_index_of_background_parameters(iparams_bg, icuts));
        end

        function fix_scale(obj)
            obj.fix_parameters(numel(obj.params))
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
                obj.set_bounds(iparams(ibnd), lb, ub, ibnd);
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
                obj.set_bounds(iparams, lb, ub, ibnd);
            end
        end

        function set_scale_bounds(obj, lb, ub)
            obj.set_bounds(size(obj.bounds,1), lb, ub);
        end

        function add_data(obj, data)
            if ~isa(data, "struct")
                data = arrayfun(@obj.convert_horace_to_struct, data);
            end
            obj.clear_cache();
            if numel(data) == 1 && isa(data.x, "cell")
                % 2D
                assert(all(isfield(data, {'x', 'y', 'e'})), ...
                       'sw_fitpowder:invalidinput', ...
                       'Input cell does not have correct fields');
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
                    assert(all(isfield(cut, {'x', 'y', 'e', 'qmin','qmax'})), ...
                           'sw_fitpowder:invalidinput', ...
                           'Input cell does not have correct fields');
                    obj.y = [obj.y cut.y(:)];
                    obj.e = [obj.e cut.e(:)];
                    dQ = (cut.qmax - cut.qmin)/obj.nQ;
                    obj.modQ_cens = [obj.modQ_cens, (cut.qmin+dQ/2):dQ:cut.qmax];
                end
            end
            obj.powspec_args.Evect = obj.ebin_cens(:)';  % row vect
        end

        function crop_energy_range(obj, emin, emax)
            % crop data
            ikeep = obj.ebin_cens >= emin & obj.ebin_cens <= emax;
            obj.apply_energy_mask(ikeep);
        end

        function exclude_energy_range(obj, elo, ehi)
            ikeep = obj.ebin_cens < elo | obj.ebin_cens > ehi;
            obj.apply_energy_mask(ikeep);
        end

        function crop_q_range(obj, qmin, qmax)
            assert(obj.ndim==2, 'sw_fitpowder:invalidinput', ...
                  '|Q| range can only be cropped on 2D datasets');
            % crop data
            ikeep = obj.modQ_cens >= qmin & obj.modQ_cens <= qmax;
            obj.modQ_cens = obj.modQ_cens(ikeep);
            obj.y = obj.y(:, ikeep);
            obj.e = obj.e(:, ikeep);
            obj.clear_cache();
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
            if obj.liveplot_interval > 0
                obj.liveplot_counter = obj.liveplot_counter + 1;
                if obj.liveplot_counter == obj.liveplot_interval
                    obj.plot_1d_or_2d(ycalc);
                    drawnow;
                    obj.liveplot_counter = 0;
                end
            end
            resid = (obj.y - ycalc);
            if obj.cost_function == "chisq"
                resid = resid./(obj.e);
            end
            % exclude nans in both ycalc and input data
            ikeep = ~isnan(resid);
            resid_sq_sum = resid(ikeep)'*resid(ikeep);
        end

        function result = fit(obj, varargin) 
            if obj.liveplot_interval > 0
                figure("color","white");
                obj.liveplot_counter = 0;
            end
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

        function estimate_constant_background(obj)
            ysort = sort(obj.y(obj.y < mean(obj.y, 'omitnan')), 'descend');
            bg = ysort(int32(numel(ysort)/2)); % default to median
            prev_skew = inf;
            for ipt = 2:numel(ysort)
                this_mean = mean(ysort(ipt:end));
                this_skew = mean((ysort(ipt:end) - this_mean).^3)/(std(ysort(ipt:end)).^3);
                if this_skew < 0
                    bg = this_mean;
                    break
                else
                    prev_skew = this_skew;
                end
            end
            if obj.background_strategy == "planar" && obj.ndim == 1
                bg = bg/obj.nQ;
            end
            % set constant background assuming last bg parameter
            obj.set_bg_parameters(obj.nparams_bg, bg);  % set constant background
        end

        function plot_result(obj, params, varargin)
            [ycalc, ~] = obj.calc_spinwave_spec(params);
            obj.plot_1d_or_2d(ycalc, varargin{:});
        end

        function plot_1d_cuts_on_data(obj, ycalc, varargin)
            modQs = mean(reshape(obj.modQ_cens, [], obj.ncuts), 1);
            for icut = 1:obj.ncuts
                ax =  subplot(1, obj.ncuts, icut);
                hold on; box on;
                plot(ax, obj.ebin_cens, obj.y(:,icut), 'ok');
                plot(ax, obj.ebin_cens, ycalc(:,icut), '-r', varargin{:});
                xlim(ax, [obj.ebin_cens(1), obj.ebin_cens(end)]);
                ymin = min(min(obj.y(:,icut)), min(ycalc(:,icut)));
                ymax = max(max(obj.y(:,icut)), max(ycalc(:,icut)));
                ylim(ax, [ymin, ymax]);
                xlabel(ax, 'Energy (meV)')
                ylabel(ax, 'Intensity');
                title(ax, num2str(modQs(icut), 2) + " $\AA^{-1}$", 'interpreter','latex')
            end
        end

        function plot_2d_contour_on_data(obj, ycalc, varargin)
            ax = subplot(1,1,1);
            box on; hold on;
            h = imagesc(ax, obj.modQ_cens, obj.ebin_cens, obj.y);
            h.AlphaData = obj.y > 0;  % make empty bins transparent
            contour(ax, obj.modQ_cens, obj.ebin_cens, ycalc, varargin{:});
            cbar = colorbar(ax);
            cbar.Label.String = "Intensity";
            xlabel(ax, "$\left|Q\right| (\AA^{-1})$", 'interpreter','latex');
            ylabel(ax, "Energy (meV)");
            ylim(ax, [obj.ebin_cens(1), obj.ebin_cens(end)]);
            xlim(ax, [obj.modQ_cens(1), obj.modQ_cens(end)]);
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
            if any(icuts == 0)
                if obj.background_strategy == "independent"
                    icuts = 1:obj.ncuts;  % apply to all cuts
                else
                    icuts = 1; % will work for planar bg
                end
            end
            % find index in vector of all parameters
            iparams = [];
            for icut = icuts
                iparams = [iparams,  iparams_bg + obj.nparams_model + (icut-1)*obj.nparams_bg];
            end
        end
        function plot_1d_or_2d(obj, ycalc, varargin)
            if obj.liveplot_interval == 0
                figure("color","white");
            else
                clf;
            end
            if obj.ndim == 1
                obj.plot_1d_cuts_on_data(ycalc, varargin{:})
            else
                obj.plot_2d_contour_on_data(ycalc, varargin{:})
            end
        end
        function apply_energy_mask(obj, ikeep)
            obj.ebin_cens = obj.ebin_cens(ikeep);
            obj.powspec_args.Evect = obj.ebin_cens(:)';
            obj.y = obj.y(ikeep, :);
            obj.e = obj.e(ikeep, :);
            obj.clear_cache();
        end
        function set_bounds(obj, iparams, lb, ub, ibnd)
            if ~isempty(lb)
                if nargin == 5
                    lb = lb(ibnd);
                end
                obj.bounds(iparams, 1) = lb;
            end
            if ~isempty(ub)
                if nargin == 5
                    ub = ub(ibnd);
                end
                obj.bounds(iparams, 2) = ub;
            end
        end
        function fix_parameters(obj, iparams)
            for iparam = iparams
               obj.bounds(iparam, :) = obj.params(iparam);
            end
        end
    end
    methods (Static=true, Hidden=true, Access = private)
        function data_struct = convert_horace_to_struct(data)
            if isa(data, "sqw")
                ndim = data.dimensions;
                assert(ndim < 3, 'spinw:fitpow:invalidinput', 'Input SQW object must be 1D or 2D');
                data_obj = data.data;
            elseif isa(data, "d1d")
                ndim = 1;
                data_obj = data;
            elseif isa(data, "d2d")
                ndim = 2;
                data_obj = data;
            else
                error('spinw:fitpow:invalidinput', 'Input must be a Horace object (sqw, d1d or d2d)')
            end
            assert(strcmp(data_obj.ulabel{1}, '|Q|'), 'spinw:fitpow:invalidinput', 'Input Horace object is not a powder cut');
            % convert edges to bin centers (need energy first if 2D)
            cens = cellfun(@(edges) (edges(1:end-1) + edges(2:end))/2, data_obj.p(ndim:-1:1), 'UniformOutput', false);
            if ndim == 2
                cens = {cens};  % otherwise MATLAB makes multiple structs 
            end
            data_struct = struct('x', cens, 'y', data_obj.s, 'e', data_obj.e);
            if ndim == 1
                % add q integration range
                data_struct.qmin = data_obj.iint(1,1);
                data_struct.qmax = data_obj.iint(2,1);
            end
        end

    end
end