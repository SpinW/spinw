classdef unittest_ndbase_estimate_hessian < sw_tests.unit_tests.unittest_super
    % Runs through unit test for @spinw/spinwave.m

    properties
        fcost = @(pars) (pars(1)+5*pars(2))^2 + (pars(1)*pars(2))^2;
        expected_hessian = @(pars) [2+2*pars(2)^2, 4*pars(1)*pars(2)+10; 
                                    4*pars(1)*pars(2)+10, 50+2*pars(1)^2];
        minimum = [0 0];
    end

    methods (Test)
        function test_wrong_number_absolute_steps(testCase)
            testCase.verifyError(...
                @() ndbase.estimate_hessian(testCase.fcost, testCase.minimum, 'step', [1,1,1]), ...
                'ndbase:estimate_hessian');
        end
        function test_automatic_step_size(testCase)
            out = ndbase.estimate_hessian(testCase.fcost, testCase.minimum);
            testCase.verify_val(out, testCase.expected_hessian(testCase.minimum), 'abs_tol', 1e-6);
        end
        function test_automatic_step_size_not_at_minimum(testCase)
            pars = [1,0];
            out = ndbase.estimate_hessian(testCase.fcost, pars);
            testCase.verify_val(out, testCase.expected_hessian(pars), 'abs_tol', 1e-6);
        end
        function test_absolute_step_size(testCase)
            pars = [1,0];
            out = ndbase.estimate_hessian(testCase.fcost, pars, 'step', pars*1e-4);
            % automatic step size assumes
            testCase.verify_val(out, testCase.expected_hessian(pars), 'abs_tol', 1e-4);
        end
        function test_relative_step_size(testCase)
            pars = [1,0];
            out = ndbase.estimate_hessian(testCase.fcost, pars, 'step', 1e-4);
            testCase.verify_val(out, testCase.expected_hessian(pars), 'abs_tol', 1e-4);
        end
        function test_absolute_step_size_negative(testCase)
            out = ndbase.estimate_hessian(testCase.fcost, testCase.minimum, 'step', [-1e-5, -1e-5]);
            testCase.verify_val(out, testCase.expected_hessian(testCase.minimum), 'abs_tol', 1e-6);
        end
        function test_warning_if_automatic_step_size_fails(testCase)
            % use function with shallow stationary point
            testCase.verifyWarning(...
                @() ndbase.estimate_hessian(@(pars) (pars(1)^3)*(pars(2)^3), ...
                                            [0,0]), ...
                'ndbase:estimate_hessian');
        end
        function test_niter(testCase)
            % set niter small so that appropriate step size not reached
            testCase.verifyWarning(...
                @() ndbase.estimate_hessian(testCase.fcost, testCase.minimum, ...
                                            'niter', 1), ...
                'ndbase:estimate_hessian');
        end
        function test_cost_tol(testCase)
            % set cost_tol so high the appropriate step size not reached
            testCase.verifyWarning(...
                @() ndbase.estimate_hessian(testCase.fcost, testCase.minimum, ...
                                            'cost_tol', 1e3), ...
                'ndbase:estimate_hessian');
        end
        function test_set_ivary(testCase)
            % use function with shallow stationary point
            fcost = @(pars) pars(1)^2 + testCase.fcost(pars(2:end));
            pars = [3, testCase.minimum];
            out = ndbase.estimate_hessian(fcost, pars, 'ivary',2:3);
            testCase.verify_val(out, testCase.expected_hessian(testCase.minimum), 'abs_tol', 1e-6);
        end
    end

end
