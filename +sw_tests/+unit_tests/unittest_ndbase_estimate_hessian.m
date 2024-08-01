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
                @() ndbase.estimate_hessian(testCase.fcost, testCase.minimum, [1,1,1]), ...
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
            out = ndbase.estimate_hessian(testCase.fcost, pars, pars*1e-4);
            % automatic step size assumes
            testCase.verify_val(out, testCase.expected_hessian(pars), 'abs_tol', 1e-6);
        end
        function test_relative_step_size(testCase)
            pars = [1,0];
            out = ndbase.estimate_hessian(testCase.fcost, pars, 1e-4);
            testCase.verify_val(out, testCase.expected_hessian(pars), 'abs_tol', 1e-6);
        end
        function test_absolute_step_size_negative(testCase)
            out = ndbase.estimate_hessian(testCase.fcost, testCase.minimum, [-1e-5, -1e-5]);
            testCase.verify_val(out, testCase.expected_hessian(testCase.minimum), 'abs_tol', 1e-6);
        end
        function test_warning_if_automatic_step_size_fails(testCase)
            % use function with shallow stationary point
            testCase.verifyWarning(...
                @() ndbase.estimate_hessian(@(pars) (pars(1)^3)*(pars(2)^3), ...
                                            [0,0]), ...
                'ndbase:estimate_hessian');
        end
    end

end
