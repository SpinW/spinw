classdef unittest_sw_fitpowder < sw_tests.unit_tests.unittest_super

    properties
        swobj = [];
        data_1d_cuts = arrayfun(@(qmin) struct('x', 1:3, 'y', 1:3, ...
                                               'e', 1:3, 'qmin', qmin, ...
                                               'qmax',qmin+1), 3.5:4.5);
        data_2d = struct('x', {{1:3, 4:5}}, 'y', [1:3; 1:3], ...
                         'e', [1:3; 1:3]);
        fit_func =  @(obj, p) matparser(obj, 'param', p, 'mat', {'J_1'}, 'init', true);
        j1 = 2;
        default_fitpow = [];
        default_fields = [];
        default_modQ_cens_1d = 3.55:0.1:5.45; % integrate over nQ pts
    end

    properties (TestParameter)
        fit_params = {{}, {'resid_handle', true}};
    end

    methods (TestClassSetup)
        function setup_spinw_obj_and_expected_result(testCase)
            % setup spinw object
            testCase.swobj = sw_model('triAF', 1);
            % subset of default fitpow object fields
            testCase.default_fitpow = struct('y', testCase.data_2d.y', ...
                                             'e', testCase.data_2d.e', ...
                                             'ebin_cens', testCase.data_2d.x{1}, ...
                                             'modQ_cens', testCase.data_2d.x{2}, ...
                                             'params', [2;0;0;0;1], ...
                                             'bounds', [-Inf   Inf;
                                                        -Inf   Inf;
                                                        -Inf   Inf;
                                                        -Inf   Inf;
                                                           0   Inf], ...
                                              'ibg', [], ...
                                              'npoly_modQ', 1, ...
                                              'npoly_en', 1);
            testCase.default_fields = fieldnames(testCase.default_fitpow);
        end
    end

    methods
        function verify_results(testCase, observed, expected, fieldnames, varargin)
            if nargin < 4
                fieldnames = testCase.default_fields;
            end
            for ifld = 1:numel(fieldnames)
                fld = fieldnames{ifld};
                testCase.verify_val(observed.(fld), expected.(fld), varargin{:});
            end
        end
    end

    methods (Test)
        function test_init_data_2d(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_2d, ...
                               testCase.fit_func, testCase.j1);
            testCase.verify_results(out, testCase.default_fitpow);
        end

        function test_init_data_1d_planar_bg(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_1d_cuts, ...
                               testCase.fit_func, testCase.j1);
            expected_fitpow = testCase.default_fitpow;
            expected_fitpow.modQ_cens = testCase.default_modQ_cens_1d;
            testCase.verify_results(out, expected_fitpow);
        end

        function test_init_data_1d_indep_bg(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_1d_cuts, ...
                               testCase.fit_func, testCase.j1, "independent");
            expected_fitpow = testCase.default_fitpow;
            expected_fitpow.modQ_cens = testCase.default_modQ_cens_1d;
            % add extra background param
            expected_fitpow.params = expected_fitpow.params([1:2,2:end],:);
            expected_fitpow.bounds = expected_fitpow.bounds([1:2,2:end],:);
            testCase.verify_results(out, expected_fitpow);
        end

        function test_set_background_strategy_to_planar(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_1d_cuts, ...
                               testCase.fit_func, testCase.j1, "independent");
            % set some background parameters for 1D cuts equivalent to a planar bg
            % with slope_en=1, slope_q=2, intercept = 3
            out.set_bg_parameters(1, 1); % en_slope = 1
            out.set_bg_parameters(2, 11, 1); % intercept = 4*2 + 3 for cut 1 
            out.set_bg_parameters(2, 13, 2); % intercept = 5*2 + 3 for cut 2
            out.set_background_strategy("planar");
            expected_fitpow = testCase.default_fitpow;
            expected_fitpow.params(2:end-1) = [2,1,3];
            expected_fitpow.modQ_cens = testCase.default_modQ_cens_1d;
            testCase.verify_results(out, expected_fitpow, ...
                                    testCase.default_fields, ...
                                    'abs_tol', 1e-10);
        end

        function test_set_background_strategy_to_indep(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_1d_cuts, ...
                               testCase.fit_func, testCase.j1, "planar");
            % set some background parameters (correpsonding to planar bg
            % with slope_en=1, slope_q=2, intercept = 3
            out.set_bg_parameters(1:3, [2,1,3]); % en_slope = 2
            out.set_background_strategy("independent");
            expected_fitpow = testCase.default_fitpow;
            % add extra background param
            expected_fitpow.params = expected_fitpow.params([1:2,2:end],:);
            expected_fitpow.bounds = expected_fitpow.bounds([1:2,2:end],:);
            expected_fitpow.params(2:2:end-1) = 1; % en slope
            expected_fitpow.params(3) = 11; % intercept = 4*2 + 3 for cut 1 
            expected_fitpow.params(5) = 13; % intercept = 5*2 + 3 for cut 2
            expected_fitpow.modQ_cens = testCase.default_modQ_cens_1d;
            testCase.verify_results(out, expected_fitpow, ...
                                    testCase.default_fields, ...
                                    'abs_tol', 1e-10);
        end

        function test_replace_2D_data_with_1D_cuts(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_2d, ...
                               testCase.fit_func, testCase.j1);
            qcens = [4, 5];
            out.replace_2D_data_with_1D_cuts(qcens-0.5, qcens+0.5)
            expected_fitpow = testCase.default_fitpow;
            expected_fitpow.modQ_cens = qcens;
            testCase.verify_results(out, expected_fitpow);
        end

        function test_replace_2D_data_with_1D_cuts_specify_bg(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_2d, ...
                               testCase.fit_func, testCase.j1);
            qcens = [4, 5];
            out.replace_2D_data_with_1D_cuts(qcens-0.5, qcens+0.5,...
                                             "independent")
            expected_fitpow = testCase.default_fitpow;
            expected_fitpow.modQ_cens = qcens;  % not nQ values as using bins in 2D data
            % add extra background param
            expected_fitpow.params = expected_fitpow.params([1:2,2:end],:);
            expected_fitpow.bounds = expected_fitpow.bounds([1:2,2:end],:);
            testCase.verify_results(out, expected_fitpow);
        end

        function test_init_data_1d_specify_nQ(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_1d_cuts, ...
                               testCase.fit_func, testCase.j1, "planar", 1);
            testCase.verify_results(out, testCase.default_fitpow);
        end

        function test_set_model_parameters(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_2d, ...
                               testCase.fit_func, 5);
            out.set_model_parameters(1, testCase.j1);
            testCase.verify_results(out, testCase.default_fitpow);
        end

        function test_set_model_parameter_bounds(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_2d, ...
                               testCase.fit_func, testCase.j1);
            lb = -5;
            ub = 5;
            out.set_model_parameter_bounds(1, -5, 5);
            expected_fitpow = testCase.default_fitpow;
            expected_fitpow.bounds(1, :) = [lb, ub];
            testCase.verify_results(out, expected_fitpow);
        end

        function test_fix_model_parameters(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_2d, ...
                               testCase.fit_func, testCase.j1);
            out.fix_model_parameters(1);
            expected_fitpow = testCase.default_fitpow;
            expected_fitpow.bounds(1, :) = [testCase.j1, testCase.j1];
            testCase.verify_results(out, expected_fitpow);
        end

        function test_set_bg_parameters_planar_bg(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_2d, ...
                               testCase.fit_func, testCase.j1);
            bg_pars = [-5, 5];
            out.set_bg_parameters(1:2, bg_pars);
            expected_fitpow = testCase.default_fitpow;
            expected_fitpow.params(2:3) = bg_pars;
            testCase.verify_results(out, expected_fitpow);
        end

        function test_set_bg_parameters_with_char(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_2d, ...
                               testCase.fit_func, testCase.j1);
            out.set_bg_parameters('E0', -5);
            expected_fitpow = testCase.default_fitpow;
            expected_fitpow.params(end-1) = -5;
            testCase.verify_results(out, expected_fitpow);
        end

        function test_set_bg_parameters_with_cell_of_char(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_2d, ...
                               testCase.fit_func, testCase.j1);
            bg_pars = [-5, 5];
            out.set_bg_parameters({'Q1', 'E1'}, bg_pars);
            expected_fitpow = testCase.default_fitpow;
            expected_fitpow.params(2:3) = bg_pars;
            testCase.verify_results(out, expected_fitpow);
        end

        function test_set_bg_parameters_indep_bg_all_cuts(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_1d_cuts, ...
                               testCase.fit_func, testCase.j1, "independent", 1);
            bg_pars = [-5, 5];
            out.set_bg_parameters(1:2, bg_pars);
            expected_fitpow = testCase.default_fitpow;
            expected_fitpow.params = [expected_fitpow.params(1);
                                      bg_pars(:); bg_pars(:); 1];
            expected_fitpow.bounds = expected_fitpow.bounds([1:2,2:end],:);
            testCase.verify_results(out, expected_fitpow);
        end

        function test_set_bg_parameters_indep_bg_specify_icut(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_1d_cuts, ...
                               testCase.fit_func, testCase.j1, "independent", 1);
            bg_pars = [-5, 5];
            out.set_bg_parameters(1:2, bg_pars, 2);
            expected_fitpow = testCase.default_fitpow;
            expected_fitpow.params = [expected_fitpow.params(1:3);
                                      bg_pars(:); 1];
            expected_fitpow.bounds = expected_fitpow.bounds([1:2,2:end],:);
            testCase.verify_results(out, expected_fitpow);
        end

        function test_fix_bg_parameters_indep_bg_specify_icut(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_1d_cuts, ...
                               testCase.fit_func, testCase.j1, "independent", 1);
            out.fix_bg_parameters(1:2, 2);
            expected_fitpow = testCase.default_fitpow;
            expected_fitpow.params = expected_fitpow.params([1:2,2:end],:);
            expected_fitpow.bounds = [expected_fitpow.bounds(1:3,:);
                                      zeros(2,2);
                                      expected_fitpow.bounds(end,:)];
            testCase.verify_results(out, expected_fitpow);
        end

        function test_fix_bg_parameters_indep_bg_all_cuts(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_1d_cuts, ...
                               testCase.fit_func, testCase.j1, "independent", 1);
            out.fix_bg_parameters(1:2);
            expected_fitpow = testCase.default_fitpow;
            expected_fitpow.params = expected_fitpow.params([1:2,2:end],:);
            expected_fitpow.bounds = [expected_fitpow.bounds(1,:);
                                      zeros(4,2);
                                      expected_fitpow.bounds(end,:)];
            testCase.verify_results(out, expected_fitpow);
        end

        function test_set_bg_parameter_bounds_indep_bg_all_cuts(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_1d_cuts, ...
                               testCase.fit_func, testCase.j1, "independent", 1);
            ub = 5;
            out.set_bg_parameter_bounds(1:2, [], [ub, ub]); % lb unchanged
            expected_fitpow = testCase.default_fitpow;
            expected_fitpow.params = expected_fitpow.params([1:2,2:end]);
            expected_fitpow.bounds = expected_fitpow.bounds([1:2,2:end],:);
            expected_fitpow.bounds(2:5, 2) = ub;
            testCase.verify_results(out, expected_fitpow);
        end

        function test_set_bg_parameter_bounds_indep_bg_specify_icut(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_1d_cuts, ...
                               testCase.fit_func, testCase.j1, "independent", 1);
            ub = 5;
            out.set_bg_parameter_bounds(1:2, [], [ub, ub], 1); % lb unchanged
            expected_fitpow = testCase.default_fitpow;
            expected_fitpow.params = expected_fitpow.params([1:2,2:end]);
            expected_fitpow.bounds = expected_fitpow.bounds([1:2,2:end],:);
            expected_fitpow.bounds(2:3, 2) = ub;
            testCase.verify_results(out, expected_fitpow);
        end
        
        function test_set_scale(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_2d, ...
                               testCase.fit_func, testCase.j1);
            scale = 5;
            out.set_scale(scale);
            expected_fitpow = testCase.default_fitpow;
            expected_fitpow.params(end) = scale;
            testCase.verify_results(out, expected_fitpow);
        end

        function test_fix_scale(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_2d, ...
                               testCase.fit_func, testCase.j1);
            out.fix_scale();
            expected_fitpow = testCase.default_fitpow;
            expected_fitpow.bounds(end, :) = 1;
            testCase.verify_results(out, expected_fitpow);
        end

        function test_set_scale_bounds(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_2d, ...
                               testCase.fit_func, testCase.j1);
            lb = 0.5;
            out.set_scale_bounds(0.5, []);
            expected_fitpow = testCase.default_fitpow;
            expected_fitpow.bounds(end, 1) = lb;
            testCase.verify_results(out, expected_fitpow);
        end

        function test_crop_energy_range_1d(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_1d_cuts, ...
                               testCase.fit_func, testCase.j1);
            out.crop_energy_range(2,5);
            expected_fitpow = testCase.default_fitpow;
            expected_fitpow.modQ_cens = testCase.default_modQ_cens_1d;
            expected_fitpow.y = expected_fitpow.y(2:end,:);
            expected_fitpow.e = expected_fitpow.e(2:end,:);
            expected_fitpow.ebin_cens = expected_fitpow.ebin_cens(2:end);
            testCase.verify_results(out, expected_fitpow);
            testCase.verify_val(out.powspec_args.Evect, ...
                                expected_fitpow.ebin_cens)
        end
        function test_crop_energy_range_2d(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_2d, ...
                               testCase.fit_func, testCase.j1);
            out.crop_energy_range(2,5);
            expected_fitpow = testCase.default_fitpow;
            expected_fitpow.y(1,:) = NaN;
            expected_fitpow.e(1,:) = NaN;
            testCase.verify_results(out, expected_fitpow);
            testCase.verify_val(out.powspec_args.Evect, ...
                                expected_fitpow.ebin_cens)
        end

        function test_exclude_energy_range_1d(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_1d_cuts, ...
                               testCase.fit_func, testCase.j1);
            out.exclude_energy_range(1.5,2.5);
            expected_fitpow = testCase.default_fitpow;
            expected_fitpow.modQ_cens = testCase.default_modQ_cens_1d;
            expected_fitpow.y = expected_fitpow.y([1,end],:);
            expected_fitpow.e = expected_fitpow.e([1,end],:);
            expected_fitpow.ebin_cens = expected_fitpow.ebin_cens([1,end]);
            testCase.verify_results(out, expected_fitpow);
            testCase.verify_val(out.powspec_args.Evect, ...
                                expected_fitpow.ebin_cens)
        end

        function test_crop_q_range(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_2d, ...
                               testCase.fit_func, testCase.j1);
            out.crop_q_range(4.5,5.5);
            expected_fitpow = testCase.default_fitpow;
            expected_fitpow.y = expected_fitpow.y(:, 2:end);
            expected_fitpow.e = expected_fitpow.e(:, 2:end);
            expected_fitpow.modQ_cens = expected_fitpow.modQ_cens(2:end);
            testCase.verify_results(out, expected_fitpow);
        end

        function test_estimate_constant_background_fails(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_2d, ...
                               testCase.fit_func, testCase.j1);
            % background estimation fails as no minimum in skew
            testCase.verifyError(...
                @() out.estimate_constant_background(), ...
                'spinw:find_indices_and_mean_of_bg_bins');
            testCase.verify_results(out, testCase.default_fitpow);
        end

        function test_estimate_constant_background(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_2d, ...
                               testCase.fit_func, testCase.j1);
            out.y(1) = 10;  % higher so other bins are background
            out.estimate_constant_background();
            expected_fitpow = testCase.default_fitpow;
            expected_fitpow.y(1) = 10;
            expected_fitpow.ibg = [3;6;2;5;4];
            expected_fitpow.params(end-1) = 2.2;
            testCase.verify_results(out, expected_fitpow);
        end

        function test_fit_background(testCase, fit_params)
            out = sw_fitpowder(testCase.swobj, testCase.data_2d, ...
                               testCase.fit_func, testCase.j1);
            out.y(1) = 10;  % higher so other bins are background
            out.fix_bg_parameters(1:2); % fix slopes of background to 0
            out.set_bg_parameters(3, 1.5); % initial guess
            out.fit_background(fit_params{:});
            expected_fitpow = testCase.default_fitpow;
            expected_fitpow.y(1) = 10;
            expected_fitpow.ibg = [3;6;2;5;4];
            expected_fitpow.params(end-1) = 2.2002;
            expected_fitpow.bounds(2:3,:) = 0; % fixed bg slopes
            testCase.verify_results(out, expected_fitpow, ...
                                    testCase.default_fields, ...
                                    'abs_tol', 1e-4);
        end

        function test_fit_background_and_scale(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_2d, ...
                               testCase.fit_func, testCase.j1);
            out.fix_bg_parameters(1:2); % fix slopes of background to 0
            out.set_bg_parameters(3, 1.5); % initial guess
            out.fit_background_and_scale();
            expected_fitpow = testCase.default_fitpow;
            expected_fitpow.params(end-1) = 0.0029;
            expected_fitpow.params(end) = 15.47;
            expected_fitpow.bounds(2:3,:) = 0; % fixed bg slopes
            testCase.verify_results(out, expected_fitpow, ...
                                    testCase.default_fields, ...
                                    'abs_tol', 1e-3);
        end

        function test_calc_cost_func_of_background_indep(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_1d_cuts, ...
                               testCase.fit_func, testCase.j1, "independent", 2);
            out.y(1) = 10;  % higher so other bins are background
            bg_pars = out.params(2:end-1);
            bg_pars(2:2:end) = 1;
            out.estimate_constant_background(); % so ibg is set
            cost = out.calc_cost_func_of_background(bg_pars);
            testCase.verify_val(cost, 10);
        end

        function test_calc_cost_func_of_background_planar(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_1d_cuts, ...
                               testCase.fit_func, testCase.j1, "planar", 2);
            out.y(1) = 10;  % higher so other bins are background
            bg_pars = out.params(2:end-1);
            bg_pars(end) = 1;
            out.estimate_constant_background(); % so ibg is set
            cost = out.calc_cost_func_of_background(bg_pars);
            testCase.verify_val(cost, 10);
        end

        function test_set_errors_of_bg_bins(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_2d, ...
                               testCase.fit_func, testCase.j1);
            out.y(1) = 10;  % high so other bins are background
            out.set_errors_of_bg_bins(100);
            expected_fitpow = testCase.default_fitpow;
            expected_fitpow.y(1) = 10;
            expected_fitpow.ibg = [3;6;2;5;4];
            expected_fitpow.e(expected_fitpow.ibg) = 100;
            testCase.verify_results(out, expected_fitpow);
        end

        function test_reset_errors_of_bg_bins(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_2d, ...
                               testCase.fit_func, testCase.j1);
            out.y(1) = 10;  % high so other bins are background
            out.set_errors_of_bg_bins(100);
            out.reset_errors_of_bg_bins();
            expected_fitpow = testCase.default_fitpow;
            expected_fitpow.y(1) = 10;
            expected_fitpow.ibg = [3;6;2;5;4];
            testCase.verify_results(out, expected_fitpow);
        end

        function test_calc_uncertainty(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_2d, ...
                               testCase.fit_func, testCase.j1);
            % perform a fit of const background only
            out.set_scale(0);
            out.powspec_args.hermit = true;
            out.fix_model_parameters(1);
            out.fix_bg_parameters(1:2);
            out.fix_scale()
            out.set_bg_parameters(3, mean(out.y, 'all'));
            [param_errors, cov] = out.calc_uncertainty(out.params);
            % test errors determined for only 1 parameter (const bg)
            testCase.verify_val(param_errors, [0; 0; 0; 0.3651; 0], 'abs_tol', 1e-4);
        end

        function test_estimate_scale_factor(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_2d, ...
                               testCase.fit_func, testCase.j1);
            out.set_bg_parameters(3, 0.05);  % set constant bg
            out.powspec_args.dE = 0.1;  % constant energy resolution
            out.powspec_args.hermit = true;
            out.estimate_scale_factor()
            expected_fitpow = testCase.default_fitpow;
            expected_fitpow.params(end) = 17.36;
            testCase.verify_results(out, expected_fitpow, ...
                                    testCase.default_fields, 'abs_tol', 1e-1);
        end

        function test_calc_cost_func_1d_planar_bg(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_1d_cuts, ...
                               testCase.fit_func, testCase.j1, "planar", 1);
            out.powspec_args.dE = 0.1;  % constant energy resolution
            out.powspec_args.hermit = true;
            cost = out.calc_cost_func(out.params);
            testCase.verify_val(cost, 25.3, 'abs_tol', 0.1);
        end

        function test_calc_cost_func_1d_indep_bg(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_1d_cuts, ...
                               testCase.fit_func, testCase.j1, "independent", 1);
            out.powspec_args.dE = 0.1;  % constant energy resolution
            out.powspec_args.hermit = true;
            cost = out.calc_cost_func(out.params);
            testCase.verify_val(cost, 25.3, 'abs_tol', 0.1);
        end

        function test_calc_cost_func_2d(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_2d, ...
                               testCase.fit_func, testCase.j1);
            out.powspec_args.dE = 0.1;  % constant energy resolution
            out.powspec_args.hermit = true;
            cost = out.calc_cost_func(out.params);
            testCase.verify_val(cost, 25.3, 'abs_tol', 0.1);
        end

        function test_add_1Dcuts_after_2D_data(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_2d, ...
                               testCase.fit_func, testCase.j1);
            out.add_data(testCase.data_1d_cuts);
            expected_fitpow = testCase.default_fitpow;
            expected_fitpow.modQ_cens = testCase.default_modQ_cens_1d; % integrtate over nQ pts
            testCase.verify_results(out, expected_fitpow);
        end

        function test_set_bg_region_data_2d(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_2d, ...
                               testCase.fit_func, testCase.j1);
            out.set_bg_region(0,1.5); % for all Q
            out.set_bg_region(2.5,3.5,4.5,inf); % for highest Q
            expected_fitpow = testCase.default_fitpow;
            expected_fitpow.ibg = [1;4;6];
            testCase.verify_results(out, expected_fitpow);
        end

        function test_set_bg_region_data_1d(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_1d_cuts, ...
                               testCase.fit_func, testCase.j1);
            expected_fitpow = testCase.default_fitpow;
            expected_fitpow.modQ_cens = testCase.default_modQ_cens_1d;
            out.set_bg_region(0,1.5); % for all cuts
            out.set_bg_region(2.5,3.5,2); % for last cut
            expected_fitpow.ibg = [1;4;6];
            testCase.verify_results(out, expected_fitpow);
        end

        function test_set_npoly_modQ(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_2d, ...
                               testCase.fit_func, testCase.j1);
            out.set_bg_npoly_modQ(2);
            expected_fitpow = testCase.default_fitpow;
            expected_fitpow.npoly_modQ = 2;
            expected_fitpow.params = [expected_fitpow.params(1:end-1); 
                                      0; 
                                      expected_fitpow.params(end)];
            expected_fitpow.bounds = [expected_fitpow.bounds(1:end-1,:); 
                                      -Inf, Inf; 
                                      expected_fitpow.bounds(end,:)];
            testCase.verify_results(out, expected_fitpow);
        end

        function test_set_npoly_en(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_2d, ...
                               testCase.fit_func, testCase.j1);
            out.set_bg_npoly_en(2);
            expected_fitpow = testCase.default_fitpow;
            expected_fitpow.npoly_en = 2;
            % add extra row in params and bounds
            expected_fitpow.params = expected_fitpow.params([1:2,2:end]); 
            expected_fitpow.bounds = expected_fitpow.bounds([1:2,2:end],:); 
            testCase.verify_results(out, expected_fitpow);
        end

        function test_set_npoly_modQ_1d_data_indep_bg(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_1d_cuts, ...
                               testCase.fit_func, testCase.j1, "independent");
            testCase.verifyError(...
                @() out.set_bg_npoly_modQ(2), ...
                'sw_fitpowder:invalidinput');
        end

        function test_set_npoly_en_1d_data_indep_bg(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_1d_cuts, ...
                               testCase.fit_func, testCase.j1, "independent");
            out.set_bg_npoly_en(2)
            expected_fitpow = testCase.default_fitpow;
            expected_fitpow.npoly_en = 2;
            expected_fitpow.modQ_cens = testCase.default_modQ_cens_1d;
            % add 3 extra rows:
            % 3 params x 2 cuts  vs 3 planar parmas order=1
            expected_fitpow.params = expected_fitpow.params([1:4,2:end]);
            expected_fitpow.bounds = expected_fitpow.bounds([1:4,2:end],:);
            testCase.verify_results(out, expected_fitpow);
        end

        function test_set_npoly_modQ_1d_data_planar_bg(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_1d_cuts, ...
                               testCase.fit_func, testCase.j1, "planar");
            % try adding order larger than numebr cuts
            testCase.verifyError(...
                @() out.set_bg_npoly_modQ(2), ...
                'sw_fitpowder:invalidinput');
            % set it to constant in modQ
            out.set_bg_npoly_modQ(0)
            expected_fitpow = testCase.default_fitpow;
            expected_fitpow.npoly_modQ = 0;
            expected_fitpow.modQ_cens = testCase.default_modQ_cens_1d;
            % remove a row
            expected_fitpow.params = expected_fitpow.params([1,3:end]);
            expected_fitpow.bounds = expected_fitpow.bounds([1,3:end],:);
            testCase.verify_results(out, expected_fitpow);
        end

        function test_export_data_2d(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_2d, ...
                               testCase.fit_func, testCase.j1);
            data = out.export_data();
            testCase.verify_val(data, testCase.data_2d);
        end

        function test_export_data_1d(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_1d_cuts, ...
                               testCase.fit_func, testCase.j1);
            data = out.export_data();
            testCase.verify_val(data, testCase.data_1d_cuts);
        end

        function test_set_bg_param_with_name_planar_bg(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_2d, ...
                               testCase.fit_func, testCase.j1);
            out.set_bg_parameters("Q1", -5);
            out.set_bg_parameters("E1", 5);
            expected_fitpow = testCase.default_fitpow;
            expected_fitpow.params(2:4) = [-5, 5,0]; % Q1, E1, E0
            testCase.verify_results(out, expected_fitpow);
        end

        function test_set_bg_param_with_invalid_name(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_2d, ...
                               testCase.fit_func, testCase.j1);
            testCase.verifyError(...
                @() out.set_bg_parameters("Q3", -5), ...
                'sw_fitpowder:invalidinput');
        end

        function test_set_bg_param_with_name_indep_bg_all_cuts(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_1d_cuts, ...
                               testCase.fit_func, testCase.j1, "independent", 1);
            bg_pars = [-5, 5];
            out.set_bg_parameters(["E1", "E0"], bg_pars);
            expected_fitpow = testCase.default_fitpow;
            expected_fitpow.params = [expected_fitpow.params(1);
                                      bg_pars(:); bg_pars(:); 1];
            expected_fitpow.bounds = expected_fitpow.bounds([1:2,2:end],:);
            testCase.verify_results(out, expected_fitpow);
        end

        function test_fix_bg_param_with_array_of_names(testCase)
            out = sw_fitpowder(testCase.swobj, testCase.data_2d, ...
                               testCase.fit_func, testCase.j1);
            bg_pars = 1:3;
            bg_labels = ["Q1","E1","E0"];
            out.set_bg_parameters(bg_labels, bg_pars);
            out.fix_bg_parameters(bg_labels);
            expected_fitpow = testCase.default_fitpow;
            expected_fitpow.params(2:4) = bg_pars;
            expected_fitpow.bounds(2:4, :) = repmat(bg_pars(:), 1,2);
            testCase.verify_results(out, expected_fitpow);
        end

    end

end
