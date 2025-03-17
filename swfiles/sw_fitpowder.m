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
% >> [pfit,cost_val,stat] = fitpow.fit('MaxIter', 1);  % passes varargin to optimizer
% >>
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
       background_strategy = "planar" % "planar" or "independent" (1D only - fbg = @(en, p)
       fbg_planar = @(en, modQ, p, npoly_en) polyval([reshape(p(1:npoly_en),1,[]), 0],en(:)) + ...
                    polyval(p(npoly_en+1:end), modQ(:)'); % en^n,...en^1, q^n,..q^1, const
       fbg_indep = @(en, p) polyval(p, en(:));
       do_cache = true
       ycalc_cached = []
       model_params_cached = []
       liveplot_counter = 0
       ibg = []
       bg_errors = []
       npoly_modQ = 1;
       npoly_en = 1;
    end

   properties (Constant)
      zero_abs_tol = 10*eps
   end
        
    methods       
        function obj = sw_fitpowder(swobj, data, fit_func, model_params, background_strategy, nQ)
            % constructor
            obj.swobj = swobj;
            obj.fit_func = fit_func;
            if nargin < 5
                background_strategy = "planar";
            elseif nargin==6
                obj.nQ = nQ; % set this before add data
            end
            obj.add_data(data)
            obj.set_background_strategy(background_strategy);
            obj.initialise_parameters_and_bounds(model_params)
        end

        function initialise_parameters_and_bounds(obj, model_params)
            obj.nparams_model = numel(model_params);
            % last parameter is scale factor (default to 1)
            obj.params = [model_params(:); 1];
            obj.bounds = [-inf, inf].*ones(numel(obj.params), 1);
            obj.bounds(end,1) = 0; % set lower bound of scale to be 0
            obj.initialise_background_parameters_and_bounds();
        end

        function initialise_background_parameters_and_bounds(obj)
            % zero intialise background parameters
            if obj.background_strategy == "independent"
                nparams_bg_total = obj.ncuts * obj.nparams_bg;
            else
                nparams_bg_total = obj.nparams_bg;
            end
            bg_params = zeros(nparams_bg_total, 1);
            bg_bounds = [-inf, inf].*ones(nparams_bg_total, 1);
            % insert in middle of params and bound
            obj.params = [obj.params(1:obj.nparams_model); bg_params; obj.params(end)];
            obj.bounds = [obj.bounds(1:obj.nparams_model,:); bg_bounds; obj.bounds(end, :)];
        end

        function set_bg_npoly_modQ(obj, npoly_modQ)
            if obj.background_strategy == "independent"
                error('sw_fitpowder:invalidinput', ...
                      'set_bg_npoly_modQ only applies to planar background');
            end
            if npoly_modQ < 0
                error('sw_fitpowder:invalidinput', ...
                      'npoly_modQ must be positive int');
            end
            if obj.ndim==1 && npoly_modQ > obj.ncuts - 1
                error('sw_fitpowder:invalidinput', ...
                      'Not enough 1D cuts to fit polynomial of order npoly_modQ.');
            end
            obj.npoly_modQ = round(npoly_modQ);
            obj.nparams_bg = obj.get_nparams_in_background_func();
            obj.initialise_background_parameters_and_bounds();
        end

        function set_bg_npoly_en(obj, npoly_en)
            if npoly_en < 0
                error('sw_fitpowder:invalidinput', ...
                      'npoly_en must be positive int');
            end
            obj.npoly_en = round(npoly_en);
            obj.nparams_bg = obj.get_nparams_in_background_func();
            obj.initialise_background_parameters_and_bounds();
        end

        function set_background_strategy(obj, strategy)
            if strategy == "planar"
                obj.fbg = @(en, modQ, p) obj.fbg_planar(en, modQ, p, obj.npoly_en);
            elseif strategy == "independent"
                if obj.ndim == 2
                    error('sw_fitpowder:invalidinput', ...
                          ['Can only have independent background strategy' ...
                          ' for 1D datasets.'])
                end
                obj.fbg = obj.fbg_indep;
            end
            % store previous bg parameters before reset size of param array
            ibg_par = (obj.nparams_model+1):(numel(obj.params)-1);
            bg_pars = obj.params(ibg_par);
            % reset param array d bounds etc.
            obj.background_strategy = strategy;
            obj.nparams_bg = obj.get_nparams_in_background_func();
            if ~isempty(obj.params)
                obj.initialise_background_parameters_and_bounds();
                if any(abs(bg_pars) > obj.zero_abs_tol) && obj.ndim == 1
                    % set good guess for background pars
                    modQs = obj.get_modQ_cens_of_cuts();
                    if obj.background_strategy == "independent"
                        % keep same energy bg
                        obj.set_bg_parameters(1:obj.npoly_en, bg_pars(1:obj.npoly_en));
                        % set intercept (zero energy) for reach cut separately
                        intercepts = obj.fbg_planar(0, modQs, bg_pars, obj.npoly_en);
                        for icut = 1:obj.ncuts
                            obj.set_bg_parameters(obj.npoly_en+1, intercepts(icut), icut);
                        end
                    elseif obj.background_strategy == "planar"
                        % reshape each col is the params for a cut
                        bg_pars = reshape(bg_pars, [], obj.ncuts);
                        % fit constant/intercepts to poly in Q
                        pval = polyfit(modQs, bg_pars(obj.npoly_en+1,:), obj.npoly_modQ);
                        obj.set_bg_parameters(obj.npoly_en+1:obj.nparams_bg, pval);
                        % average energy dep. params (only very rough)
                        bg_pars_en = mean(bg_pars(1:obj.npoly_en,:), 2);
                        obj.set_bg_parameters(1:obj.npoly_en, bg_pars_en);
                    end
                end
            end
        end

        function fix_model_parameters(obj, iparams)
            if any(iparams > obj.nparams_model)
                error('sw_fitpowder', 'Parameter indices supplied must be within number of model parameters');
            end
            obj.fix_parameters(iparams)
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
            obj.params(iparams) = values;
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
                if obj.ndim == 2
                    % clear previously added 2D data
                    obj.modQ_cens = [];
                    obj.y = [];
                    obj.e = [];
                    warning('spinw:sw_fitpowder:add_data', ...
                            'Clearing 2D data previously added.');
                end
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

        function replace_2D_data_with_1D_cuts(obj, qmins, qmaxs, background_strategy)
            assert(numel(qmins)==numel(qmaxs), 'sw_fitpowder:invalidinput', ...
                  'Must pass same number of mins and maxs to make cuts.');
            cuts = [];
            for icut = 1:numel(qmins)
                cuts = [cuts obj.cut_2d_data(qmins(icut), qmaxs(icut))];
            end
            obj.add_data(cuts);
            if nargin > 3 && background_strategy ~= obj.background_strategy
                obj.set_background_strategy(background_strategy)
            end
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

        function bg = calc_background(obj, bg_params)
            % add background
            if obj.background_strategy == "planar"
                bg = obj.fbg(obj.ebin_cens, obj.modQ_cens, bg_params);
            elseif obj.ndim == 1
                % add energy dependent background to each cut
                bg = zeros(size(obj.y));
                istart = 1;
                for icut = 1:obj.ncuts
                    iend = istart + obj.nparams_bg - 1;
                    bg(:,icut) = obj.fbg(obj.ebin_cens(:), bg_params(istart:iend));
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
            obj.clear_background_region();
        end

        function clear_background_region(obj)
            obj.ibg = [];
            obj.reset_errors_of_bg_bins();
        end

        function set_bg_region(obj, en_lo, en_hi, varargin)
            % 2D data: obj.set_bg_region(en_lo, en_hi, q_lo, q_hi)
            % 1D data: obj.set_bg_region(en_lo, en_hi, icuts)
            % both:
            % obj.set_bg_region(en_lo, en_hi)  % apply to all cuts/Q
            if isempty(varargin)
                iq = 1:size(obj.y, 2);
            elseif obj.ndim == 1 && numel(varargin)==1
                iq = varargin{1};
            elseif obj.ndim==2 && numel(varargin)==2
                [q_lo, q_hi] = varargin{:};
                iq = obj.modQ_cens < q_hi & obj.modQ_cens > q_lo;
            else
                error('spinw:sw_fitpowder:set_bg_region', ...
                      ['Wrong number of additional aruguments: 2 required' ...
                      ' for 2D problem (qlo and qhi)' ...
                      'and 1 required for 1D problem (icut indices)']);
            end
            ien = obj.ebin_cens < en_hi & obj.ebin_cens > en_lo;
            mask = false(size(obj.y));
            mask(ien, iq) = true;
            mask = mask & isfinite(obj.y);
            obj.ibg = [obj.ibg; find(mask(:))];
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
            bg = calc_background(obj, params(obj.nparams_model + 1:end-1));
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

        function resid = calc_residuals(obj, params)
            % evaluate fit function
            [ycalc, ~] = obj.calc_spinwave_spec(params);
            if obj.liveplot_interval > 0
                obj.liveplot_counter = obj.liveplot_counter + 1;
                if obj.liveplot_counter == obj.liveplot_interval
                    obj.plot_1d_or_2d(ycalc);
                    drawnow;
                    obj.liveplot_counter = 0;
                end
            end
            resid = obj.eval_residuals(obj.y, obj.e, ycalc);
        end

        function resid_sq_sum = calc_cost_func(obj, params)
            resid = obj.calc_residuals(params);
            resid_sq_sum = resid'*resid;
        end

        function resid = calc_residuals_of_background(obj, bg_params)
            bg = obj.calc_background(bg_params);
            if obj.ndim == 1 && obj.background_strategy=="planar"
                % integrate nQ |Q| points for each cut
                bg = obj.rebin_powspec_to_1D_cuts(bg);
            end
            resid = obj.eval_residuals(obj.y(obj.ibg), obj.e(obj.ibg), bg(obj.ibg));
        end

        function resid_sq_sum = calc_cost_func_of_background(obj, bg_params)
            resid = obj.calc_residuals_of_background(bg_params);
            resid_sq_sum = resid'*resid;
        end

        function varargout = fit(obj, varargin) 
            if obj.liveplot_interval > 0
                figure("color","white");
                obj.liveplot_counter = 0;
            end
            % setup cell for output of ndbase optimizer/minimizer
            varargout = cell(1,nargout(obj.optimizer));
            % pass params and bounds as rows for ndbase optimisers
            if any(cellfun(@(elem) elem =="resid_handle", varargin(1:2:end)))
                fobj = @obj.calc_residuals; % minimise using residual array
            else
                fobj = @obj.calc_cost_func; % minimise scalar
            end
            [varargout{:}] = obj.optimizer([], fobj, obj.params(:)', ...
                                       'lb', obj.bounds(:,1)', ...
                                       'ub', obj.bounds(:,2)', ...
                                       varargin{:});
        end

        function varargout = fit_background_and_scale(obj, varargin) 
            % fix all model parameters
            initial_bounds = obj.bounds; % store so bounds can be reset
            obj.fix_model_parameters(1:obj.nparams_model);
            varargout = cell(1,nargout(obj.optimizer));
            [varargout{:}] = obj.fit(varargin{:});
            % reset bounds
            obj.bounds = initial_bounds;
            % overwrite non model parameters
            obj.params = varargout{1}(:); % assume first output (as for ndbase)
        end


        function [param_errors, varargout] = calc_uncertainty(obj, params, varargin)
        % Function to estimate parameter errors by evaluating hessian.
        % By default will only evaluate derivatives for parameters
        % which are free to vary in the fit. The index of the
        % parameters is optionally output. The second optional output
        % is the full covariance matrix.
        %
        % ### Usage
        %
        % errors = sw_fitpowder_obj.calc_uncertainty(params)
        % [errors, cov] = sw_fitpowder_obj.calc_uncertainty(params)
        %
        % ### See Also
        % 
        % [ndbase.estimate_hessian]
        %
            if ~any(cellfun(@(elem) elem =="ivary", varargin(1:2:end)))
                % set parameters varied in fit from bounds
                ivary = find(obj.bounds(:,2) - obj.bounds(:,1) > obj.zero_abs_tol);
                varargin(end+1:end+2) = {"ivary", ivary};
            end
            [hess, stats] = ndbase.estimate_hessian(@obj.calc_cost_func, params, varargin{:});
            % determine num DoF
            nbins = sum(isfinite(obj.y) & obj.e > obj.zero_abs_tol, 'all');
            ndof = nbins - size(hess,1);
            % calc errors
            cov = inv(hess) * 2.0 * stats.cost_val/ ndof;
            param_errors = zeros(size(params));
            param_errors(ivary) = sqrt(diag(cov));
            varargout = {cov}; 
        end

        function varargout = fit_background(obj, varargin)
            if isempty(obj.ibg)
                obj.find_indices_and_mean_of_bg_bins();
            end
            % check if enough bins for parameters
            ibg_par = (obj.nparams_model+1):(numel(obj.params)-1);
            lb = obj.bounds(ibg_par,1)';
            ub = obj.bounds(ibg_par,2)';
            ifixed = sum(abs(lb - ub) < obj.zero_abs_tol);
            if numel(obj.ibg) < numel(ibg_par) - ifixed
                error('spinw:sw_fitpowder:fit_background', ...
                      'Not enough points to fit the function.');
            end
            if any(cellfun(@(elem) elem =="resid_handle", varargin(1:2:end)))
                fobj = @obj.calc_residuals_of_background; % minimise using residual array
            else
                fobj = @obj.calc_cost_func_of_background; % minimise scalar
            end
            varargout = cell(1,nargout(obj.optimizer));
            [varargout{:}] = obj.optimizer([], fobj, ...
                                           obj.params(ibg_par)', ...
                                          'lb', obj.bounds(ibg_par,1)', ...
                                          'ub', obj.bounds(ibg_par,2)', ...
                                           varargin{:});
            obj.params(ibg_par) = varargout{1}(:);
        end

        function estimate_scale_factor(obj)
            % set scale factor to 1
            params = obj.params;
            params(end) = 1;
            [ycalc, bg] = calc_spinwave_spec(obj, params);
            if obj.ndim == 1 && obj.background_strategy=="planar"
                % integrate nQ |Q| points for each cut
                bg = obj.rebin_powspec_to_1D_cuts(bg);
            end
            scale = max(obj.y - bg, [], "all")/max(ycalc - bg, [], "all");
            obj.params(end) = scale;
        end

        function estimate_constant_background(obj)
            if isempty(obj.ibg)
                bg = obj.find_indices_and_mean_of_bg_bins();
            else
                bg = mean(obj.y(obj.ibg));
            end
            % set constant background assuming last bg parameter
            obj.set_bg_parameters(obj.nparams_bg, bg);  % set constant background
        end

        function set_errors_of_bg_bins(obj, val)
            if isempty(obj.ibg)
                obj.find_indices_and_mean_of_bg_bins();
            end
            if nargin < 2
                val = max(obj.e(obj.ibg));
            end
            obj.bg_errors = obj.e(obj.ibg); % save values so can reset
            obj.e(obj.ibg) = val;
        end

        function reset_errors_of_bg_bins(obj)
            if ~isempty(obj.ibg) && ~isempty(obj.bg_errors)
                obj.e(obj.ibg) = obj.bg_errors;
            end
            obj.bg_errors = [];
        end

        function plot_result(obj, params, varargin)
            [ycalc, ~] = obj.calc_spinwave_spec(params);
            obj.plot_1d_or_2d(ycalc, varargin{:});
        end

        function plot_1d_cuts_on_data(obj, ycalc, varargin)
            modQs = obj.get_modQ_cens_of_cuts();
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
            ax = obj.plot_2d_data(obj.y);
            contour(ax, obj.modQ_cens, obj.ebin_cens, ycalc, varargin{:});
            legend('Calculated');
        end

        function ax = plot_2d_data(obj, y)
            ax = gca;
            box on; hold on;
            h = imagesc(ax, obj.modQ_cens, obj.ebin_cens, y);
            h.AlphaData = double(obj.y > 0);  % make empty bins transparent
            cbar = colorbar(ax);
            cbar.Label.String = "Intensity";
            xlabel(ax, "$\left|Q\right| (\AA^{-1})$", 'interpreter','latex');
            ylabel(ax, "Energy (meV)");
            ylim(ax, [obj.ebin_cens(1), obj.ebin_cens(end)]);
            xlim(ax, [obj.modQ_cens(1), obj.modQ_cens(end)]);
        end

        function plot_background_region(obj)
            if ~isempty(obj.ibg)
                if obj.ndim == 1
                    obj.plot_background_region_1d()
                else
                    obj.plot_background_region_2d()
                end
            end
        end

        function plot_1d_cuts_of_2d_data(obj, qmins, qmaxs, params)
            % optionally plot ycalc provided, otherwise will plot
            % fitpow.ycalc if not empty
            assert(obj.ndim ==2, ...
                   'sw_fitpowder:invalidinput', ...
                   'This function is only valid for 2D data');
            if nargin > 3
                [ycalc, ~] = obj.calc_spinwave_spec(params);
            else
                ycalc = [];
            end
            figure("color","white");
            ncuts = numel(qmins);
            for icut = 1:ncuts
                ikeep = obj.modQ_cens > qmins(icut) & obj.modQ_cens <= qmaxs(icut);
                ax =  subplot(1, ncuts, icut);
                hold on; box on;
                ycut = sum(obj.y(:,ikeep), 2);
                plot(ax, obj.ebin_cens, ycut, 'ok');
                if ~isempty(ycalc)
                    plot(ax, obj.ebin_cens, sum(ycalc(:,ikeep), 2), '-r');
                end
                % calc xlims
                ifinite = isfinite(ycut);
                istart = find(ifinite, 1, 'first');
                iend = find(ifinite, 1, 'last');
                xlim(ax, [obj.ebin_cens(istart), obj.ebin_cens(iend)]);
                xlabel(ax, 'Energy (meV)')
                ylabel(ax, 'Intensity');
                title(ax, num2str(0.5*(qmaxs(icut)+qmins(icut)), 2) + " $\AA^{-1}$", 'interpreter','latex')
            end
        end

    end

    % private
    methods (Hidden=true, Access = private)
        function resid = eval_residuals(obj, y, e, ycalc)
            resid = (y - ycalc);
            if obj.cost_function == "chisq"
                resid = resid./e;
            end
            % exclude nans in both ycalc and input data
            resid = resid(isfinite(resid) & e > obj.zero_abs_tol);
        end


        function modQ_cens = get_modQ_cens_of_cuts(obj)
            modQ_cens = mean(reshape(obj.modQ_cens, [], obj.ncuts), 1);
        end

        function cut = cut_2d_data(obj, qmin, qmax)
            assert(obj.ndim ==2, ...
                   'sw_fitpowder:invalidinput', ...
                   'This function is only valid for 2D data');
            cut = struct('x', obj.ebin_cens, 'qmin',  qmin, 'qmax', qmax);
            ikeep = obj.modQ_cens > qmin & obj.modQ_cens <= qmax;
            ifinite = isfinite(obj.y(:, ikeep));
            cut.y = mean(obj.y(:, ikeep), 2, 'omitnan');
            cut.e = sqrt(sum(obj.e(:, ikeep).^2, 2))./sum(ifinite, 2);
        end
        function ycalc = rebin_powspec_to_1D_cuts(obj, ycalc)
            % sum up successive nQ points along |Q| axis (dim=2)
            ycalc = reshape(ycalc, size(ycalc,1), obj.nQ, []);
            ycalc = squeeze(mean(ycalc, 2, 'omitnan'));
        end

        function nbg_pars = get_nparams_in_background_func(obj)
            % get background parameters
            if obj.background_strategy == "planar"
                nbg_pars = obj.npoly_en + obj.npoly_modQ + 1; % dependent variables energy, modQ
            else
                nbg_pars = obj.npoly_en + 1; % dependent variable energy
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
            if obj.ndim == 1
                obj.ebin_cens = obj.ebin_cens(ikeep);
                obj.powspec_args.Evect = obj.ebin_cens(:)';
                obj.y = obj.y(ikeep, :);
                obj.e = obj.e(ikeep, :);
            else
                obj.y(~ikeep, :) = NaN;
                obj.e(~ikeep, :) = NaN;
            end
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
        function bg = find_indices_and_mean_of_bg_bins(obj)
             % start with a seed of indices of likely non-signal bins
            iseed = find(isfinite(obj.y) & obj.y < mean(obj.y(:), 'omitnan'));
            [ysort, isort] = sort(obj.y(iseed), 'descend');
            bg = NaN;
            prev_skew = inf;
            for ipt = 1:numel(isort)
                this_mean = mean(ysort(ipt:end));
                this_skew = mean((ysort(ipt:end) - this_mean).^3);
                if this_skew < 0 || this_skew > prev_skew
                    bg = this_mean;
                    break
                else
                    prev_skew = this_skew;
                end
            end
            % store indices of background bins
            if ~isfinite(bg)
                error('spinw:find_indices_and_mean_of_bg_bins', ...
                      'Could not estimate background.');
            else
                obj.ibg = iseed(isort(ipt:end));
            end
        end

        function plot_background_region_2d(obj)
            im = findobj(gca, 'type', 'Image');
            is_valid = ~isempty(im) &&  all(size(im.CData) == size(obj.y));
            assert(is_valid, 'sw_fitpowder:invalidinput', ...
                   ['This function requires active axes to have an' ...
                   ' image/colorfill plot of the data.']);
            im.AlphaData(obj.ibg) = 0.5;
        end

        function plot_background_region_1d(obj)
            fig = gcf();
            axs = flip(fig.Children);
            assert(obj.ncuts==numel(axs), 'sw_fitpowder:invalidinput', ...
                   ['This function requires active figure to have same' ...
                   ' number of axes as 1D cuts']);
            [ien, icuts] = ind2sub(size(obj.y), obj.ibg);
            for icut = 1:obj.ncuts
                ibg_cut = ien(icuts == icut);
                plot(axs(icut), obj.ebin_cens(ibg_cut), ...
                                obj.y(ibg_cut,icut), 'xb')
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