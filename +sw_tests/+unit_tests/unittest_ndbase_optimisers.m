classdef unittest_ndbase_optimisers < sw_tests.unit_tests.unittest_super
    % Runs through unit test for ndbase optimisers using bounded parameter
    % transformations.

    properties
        rosenbrock = @(x) (1-x(1)).^2 + 100*(x(2) - x(1).^2).^2;
        rosenbrock_minimum = [1, 1];
    end

    properties (TestParameter)
        optimiser = {@ndbase.simplex, @ndbase.lm4};
        poly_func = {@(x, p) polyval(p, x), '@(x, p) polyval(p, x)'}
    end

    methods (Test)
        function test_optimise_data_struct(testCase, optimiser, poly_func)
            linear_pars = [2, 1];
            dat = struct('x', 1:3, 'e', ones(1,3));
            dat.y = polyval(linear_pars, dat.x);
            [pars_fit, cost_val, ~] = optimiser(dat, poly_func, [-1,-1]);
            testCase.verify_val(pars_fit, linear_pars, 'abs_tol', 1e-3);
            testCase.verify_val(cost_val, 0, 'abs_tol', 1e-6);
        end

        function test_optimise_residual_array_lm(testCase, optimiser)
            linear_pars = [2, 1];
            x = 1:3;
            y = polyval(linear_pars, x);
            [pars_fit, cost_val, ~] = optimiser([], @(p) y - polyval(p, x), [-1,-1], 'resid_handle', true);
            testCase.verify_val(pars_fit, linear_pars, 'abs_tol', 1e-3);
            testCase.verify_val(cost_val, 0, 'abs_tol', 1e-6);
        end

        function test_optimise_rosen_free(testCase, optimiser)
            [pars_fit, cost_val, ~] = optimiser([], testCase.rosenbrock, [-1,-1]);
            testCase.verify_val(pars_fit, testCase.rosenbrock_minimum, 'abs_tol', 1e-3);
            testCase.verify_val(cost_val, 0, 'abs_tol', 1e-6);
        end

        function test_optimise_rosen_lower_bound_minimum_accessible(testCase, optimiser)
            [pars_fit, cost_val, ~] = optimiser([], testCase.rosenbrock, [-1,-1], 'lb', [-2, -2]);
            testCase.verify_val(pars_fit, testCase.rosenbrock_minimum, 'abs_tol', 1e-3);
            testCase.verify_val(cost_val, 0, 'abs_tol', 1e-6);
        end

        function test_optimise_rosen_lower_bound_minimum_not_accessible(testCase, optimiser)
            % note intital guess is outside bounds
            [pars_fit, cost_val, ~] = optimiser([], testCase.rosenbrock, [-1,-1], 'lb', [-inf, 2]);
            testCase.verify_val(pars_fit, [-1.411, 2], 'abs_tol', 1e-3);
            testCase.verify_val(cost_val, 5.821, 'abs_tol', 1e-3);
        end

        function test_optimise_rosen_upper_bound_minimum_accessible(testCase, optimiser)
            [pars_fit, cost_val, ~] = optimiser([], testCase.rosenbrock, [-1,-1], 'ub', [2, 2]);
            testCase.verify_val(pars_fit, testCase.rosenbrock_minimum, 'abs_tol', 1e-3);
            testCase.verify_val(cost_val, 0, 'abs_tol', 1e-6);
        end

        function test_optimise_rosen_upper_bound_minimum_not_accessible(testCase, optimiser)
            [pars_fit, cost_val, ~] = optimiser([], testCase.rosenbrock, [-1,-1], 'ub', [0, inf]);
            testCase.verify_val(pars_fit, [0, 0], 'abs_tol', 1e-3);
            testCase.verify_val(cost_val, 1, 'abs_tol', 2e-3);
        end

        function test_optimise_rosen_both_bounds_minimum_accessible(testCase, optimiser)
            [pars_fit, cost_val, ~] = optimiser([], testCase.rosenbrock, [-1,-1], 'lb', [-5, -5], 'ub', [5, 5]);
            testCase.verify_val(pars_fit, testCase.rosenbrock_minimum, 'abs_tol', 1e-3);
            testCase.verify_val(cost_val, 0, 'abs_tol', 1e-6);
        end

        function test_optimise_rosen_both_bounds_minimum_not_accessible(testCase)
            % note intital guess is outside bounds
            [pars_fit, cost_val, ~] = ndbase.simplex([], testCase.rosenbrock, [-1,-1], 'lb', [-0.5, -0.5], 'ub', [0, 0]);
            testCase.verify_val(pars_fit, [0, 0], 'abs_tol', 1e-3);
            testCase.verify_val(cost_val, 1, 'abs_tol', 1e-6);
        end

        function test_optimise_rosen_parameter_fixed_minimum_not_accessible(testCase)
            % note intital guess is outside bounds
            [pars_fit, cost_val, ~] = ndbase.simplex([], testCase.rosenbrock, [-1,-1], 'lb', [0, -0.5], 'ub', [0, 0]);
            testCase.verify_val(pars_fit, [0, 0], 'abs_tol', 1e-3);
            testCase.verify_val(cost_val, 1, 'abs_tol', 1e-6);
        end

        function test_optimise_rosen_parameter_all_fixed(testCase, optimiser)
            % note intital guess is outside bounds
            [pars_fit, cost_val, ~] = optimiser([], testCase.rosenbrock, [-1,-1], 'lb', [0, 0], 'ub', [0, 0]);
            testCase.verify_val(pars_fit, [0, 0], 'abs_tol', 1e-3);
            testCase.verify_val(cost_val, 1, 'abs_tol', 1e-6);
        end

    end
end