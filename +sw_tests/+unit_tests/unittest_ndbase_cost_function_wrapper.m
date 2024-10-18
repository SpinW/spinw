classdef unittest_ndbase_cost_function_wrapper < sw_tests.unit_tests.unittest_super
    % Runs through unit test for ndbase optimisers, atm only simplex passes
    % these tests

    properties
        fcost = @(p) (p(1)-1)^2 + (p(2)-2)^2
        params = [2,4]
    end

    properties (TestParameter)
        bound_param_name = {'lb', 'ub'}
        no_lower_bound = {[], [-inf, -inf]};
        no_upper_bound = {[], [inf, inf]};
        errors = {ones(1,3), [], zeros(1,3), 'NoField'}
    end

    methods
        function [pfree, pbound, cost_val] = get_pars_and_cost_val(testCase, cost_func_wrap)
            pfree = cost_func_wrap.get_free_parameters(testCase.params);
            pbound = cost_func_wrap.get_bound_parameters(pfree);
            cost_val = cost_func_wrap.eval_cost_function(pfree);
        end
    end

    methods (Test)

        function test_init_with_fcost_no_bounds(testCase)
            cost_func_wrap = ndbase.cost_function_wrapper(testCase.fcost, testCase.params);
            [pfree, pbound, cost_val] = testCase.get_pars_and_cost_val(cost_func_wrap);
            testCase.verify_val(pfree, testCase.params);
            testCase.verify_val(pbound, testCase.params);
            testCase.verify_val(cost_val, testCase.fcost(pbound), 'abs_tol', 1e-4);
        end

        function test_init_with_fcost_no_bounds_name_value_passed(testCase, no_lower_bound, no_upper_bound)
            cost_func_wrap = ndbase.cost_function_wrapper(testCase.fcost, testCase.params, 'lb', no_lower_bound, 'ub', no_upper_bound);
            [pfree, pbound, cost_val] = testCase.get_pars_and_cost_val(cost_func_wrap);
            testCase.verify_val(pfree, testCase.params);
            testCase.verify_val(pbound, testCase.params);
            testCase.verify_val(cost_val, testCase.fcost(pbound), 'abs_tol', 1e-4);
        end

        function test_init_with_fcost_lower_bound_only(testCase)
            % note first param outside bounds
            cost_func_wrap = ndbase.cost_function_wrapper(testCase.fcost, testCase.params, 'lb', [3, 1]);
            [pfree, pbound, cost_val] = testCase.get_pars_and_cost_val(cost_func_wrap);
            testCase.verify_val(pfree, [0, 3.8730], 'abs_tol', 1e-4);
            testCase.verify_val(pbound, [3, 4], 'abs_tol', 1e-4);
            testCase.verify_val(cost_val, testCase.fcost(pbound), 'abs_tol', 1e-4);
        end

        function test_init_with_fcost_upper_bound_only(testCase)
            % note second param outside bounds
            cost_func_wrap = ndbase.cost_function_wrapper(testCase.fcost, testCase.params, 'ub', [3, 1]);
            [pfree, pbound, cost_val] = testCase.get_pars_and_cost_val(cost_func_wrap);
            testCase.verify_val(pfree, [1.7320, 0], 'abs_tol', 1e-4);
            testCase.verify_val(pbound, [2, 1], 'abs_tol', 1e-4);
            testCase.verify_val(cost_val, testCase.fcost(pbound), 'abs_tol', 1e-4);
        end

        function test_init_with_fcost_both_bounds(testCase)
            % note second param outside bounds
            cost_func_wrap = ndbase.cost_function_wrapper(testCase.fcost, testCase.params, 'lb', [1, 2], 'ub', [3, 2.5]);
            [pfree, pbound, cost_val] = testCase.get_pars_and_cost_val(cost_func_wrap);
            testCase.verify_val(pfree, [0, 1.5708], 'abs_tol', 1e-4);
            testCase.verify_val(pbound, [2, 2.5], 'abs_tol', 1e-4);
            testCase.verify_val(cost_val, testCase.fcost(pbound), 'abs_tol', 1e-4);
        end

        function test_init_with_fcost_both_bounds_with_fixed_param(testCase)
            % note second param outside bounds
            cost_func_wrap = ndbase.cost_function_wrapper(testCase.fcost, testCase.params, 'lb', [1, 2.5], 'ub', [3, 2.5]);
            [pfree, pbound, cost_val] = testCase.get_pars_and_cost_val(cost_func_wrap);
            testCase.verify_val(pfree, 0, 'abs_tol', 1e-4);  % only first param free
            testCase.verify_val(pbound, [2, 2.5], 'abs_tol', 1e-4);
            testCase.verify_val(cost_val, testCase.fcost(pbound), 'abs_tol', 1e-4);
            testCase.verify_val(cost_func_wrap.ifixed, 2);
            testCase.verify_val(cost_func_wrap.ifree, 1);
            testCase.verify_val(cost_func_wrap.pars_fixed, 2.5);
        end


        function test_init_with_fcost_both_bounds_with_fixed_param_using_ifix(testCase)
            % note second param outside bounds
            cost_func_wrap = ndbase.cost_function_wrapper(testCase.fcost, testCase.params, 'lb', [1, 2], 'ub', [3, 2.5], 'ifix', [2]);
            [pfree, pbound, cost_val] = testCase.get_pars_and_cost_val(cost_func_wrap);
            testCase.verify_val(pfree, 0, 'abs_tol', 1e-4);  % only first param free
            testCase.verify_val(pbound, [2, 2.5], 'abs_tol', 1e-4);
            testCase.verify_val(cost_val, testCase.fcost(pbound), 'abs_tol', 1e-4);
            testCase.verify_val(cost_func_wrap.ifixed, 2);
            testCase.verify_val(cost_func_wrap.ifree, 1);
            testCase.verify_val(cost_func_wrap.pars_fixed, 2.5);
        end

        function test_init_with_fcost_no_bounds_with_fixed_param_using_ifix(testCase)
            % note second param outside bounds
            cost_func_wrap = ndbase.cost_function_wrapper(testCase.fcost, testCase.params, 'ifix', [2]);
            [pfree, pbound, cost_val] = testCase.get_pars_and_cost_val(cost_func_wrap);
            testCase.verify_val(pfree, testCase.params(1), 'abs_tol', 1e-4);  % only first param free
            testCase.verify_val(pbound, testCase.params, 'abs_tol', 1e-4);
            testCase.verify_val(cost_func_wrap.ifixed, 2);
            testCase.verify_val(cost_func_wrap.ifree, 1);
            testCase.verify_val(cost_func_wrap.pars_fixed, testCase.params(2));
        end

        function test_init_with_data(testCase, errors)
            % all errors passed lead to unweighted residuals (either as
            % explicitly ones or the default weights if invalid errors)
            if ischar(errors) && errors == "NoField"
                dat = struct('x', 1:3);
            else
                dat = struct('x', 1:3, 'e', errors);
            end
            dat.y = polyval(testCase.params, dat.x);
            cost_func_wrap = ndbase.cost_function_wrapper(@(x, p) polyval(p, x), testCase.params, 'data', dat);
            [pfree, pbound, cost_val] = testCase.get_pars_and_cost_val(cost_func_wrap);
            testCase.verify_val(pfree, testCase.params, 'abs_tol', 1e-4);
            testCase.verify_val(pbound, testCase.params, 'abs_tol', 1e-4);
            testCase.verify_val(cost_val, 0, 'abs_tol', 1e-4);
        end

        function test_wrong_size_bounds(testCase, bound_param_name)
            testCase.verifyError(...
                @() ndbase.cost_function_wrapper(testCase.fcost, testCase.params, bound_param_name, ones(3)), ...
                'ndbase:cost_function_wrapper:WrongInput');
        end

        function test_incompatible_bounds(testCase)
            testCase.verifyError(...
                @() ndbase.cost_function_wrapper(testCase.fcost, testCase.params, 'lb', [1,1,], 'ub', [0,0]), ...
                'ndbase:cost_function_wrapper:WrongInput');
        end
        
    end
end