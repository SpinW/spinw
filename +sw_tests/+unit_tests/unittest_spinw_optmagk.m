classdef unittest_spinw_optmagk < sw_tests.unit_tests.unittest_super

    properties
        swobj = [];
        % output from optmagk
        default_mag_str = struct('nExt', int32([1 1 1]), ...
                                 'F', [sqrt(1/3) + 1i*sqrt(1/2); ...
                                       sqrt(1/3); ...
                                       sqrt(1/3) - 1i*sqrt(1/2)], ...
                                 'k', [1; 0; 0]);
        orig_rng_state = [];
    end
    properties (TestParameter)
        kbase_opts = {[1; 1; 0], [1 0; 0 1; 0 0]};
    end
    methods (TestClassSetup)
        function set_seed(testCase)
            testCase.orig_rng_state = rng;
            rng('default');
        end
    end
    methods (TestMethodSetup)
        function setup_chain_model(testCase)            
            testCase.swobj = spinw();
            testCase.swobj.genlattice('lat_const', [3 8 8])
            testCase.swobj.addatom('r',[0; 0; 0],'S',1)
            testCase.swobj.gencoupling();
        end
    end
    methods(TestMethodTeardown)
        function reset_seed(testCase)
            rng(testCase.orig_rng_state);
        end
    end
    methods (Test, TestTags = {'Symbolic'})
        function test_symbolic_warns_returns_nothing(testCase)
            testCase.swobj.addmatrix('label', 'J1', 'value', 1);
            testCase.swobj.addcoupling('mat', 'J1', 'bond', 1);
            testCase.swobj.symbolic(true)
            testCase.verifyWarning(...
                @() testCase.swobj.optmagk, ...
                'spinw:optmagk:NoSymbolic');
            testCase.verifyEmpty(testCase.swobj.mag_str.k);
            testCase.verifyEmpty(testCase.swobj.mag_str.F);
            testCase.verify_val(testCase.swobj.mag_str.nExt, ...
                                int32([1 1 1]));
        end
    end
    methods (Test)
        function test_wrong_shape_kbase_raises_error(testCase)
            testCase.verifyError(...
                @() testCase.swobj.optmagk('kbase', [1 1 0]), ...
                'spinw:optmagk:WrongInput');
        end
        function test_fm_chain_optk(testCase)
            testCase.swobj.addmatrix('label', 'J1', 'value', -1);
            testCase.swobj.addcoupling('mat', 'J1', 'bond', 1);
            out = testCase.swobj.optmagk('seed', 1);
            out.stat = rmfield(out.stat, 'nFunEvals');

            expected_mag_str = testCase.default_mag_str;
            expected_out = struct('k', expected_mag_str.k, ...
                                  'E', -1, ....
                                  'F', expected_mag_str.F, ...
                                  'stat', struct('S', 0, ...
                                                 'exitflag', -1));
            % Test struct output by optmagk
            testCase.verify_val(out, expected_out, 'abs_tol', 2e-4);
            % Also test spinw attributes have been set
            expected_mag_str = testCase.default_mag_str;
            testCase.verify_val(testCase.swobj.mag_str, expected_mag_str, ...
                                'abs_tol', 2e-4);
        end
        function test_afm_chain_optk(testCase)
            testCase.swobj.addmatrix('label', 'J1', 'value', 1);
            testCase.swobj.addcoupling('mat', 'J1', 'bond', 1);
            % Use seed for reproducibility
            testCase.swobj.optmagk('seed', 1);
            expected_mag_str = testCase.default_mag_str;
            expected_mag_str.k = [0.5; 0; 0];
            testCase.verify_val(testCase.swobj.mag_str, expected_mag_str, ...
                                'abs_tol', 1e-4);
        end
        function test_kbase(testCase, kbase_opts)
            % See https://doi.org/10.1103/PhysRevB.59.14367
            swobj = spinw();
            swobj.genlattice('lat_const', [3 3 8])
            swobj.addatom('r',[0; 0; 0],'S',1)
            swobj.gencoupling();
            J1 = 1.2;
            J2 = 1.0;
            swobj.addmatrix('label', 'J1', 'value', J1);
            swobj.addmatrix('label', 'J2', 'value', J2);
            swobj.addcoupling('mat', 'J1', 'bond', 2, 'subidx', 2);
            swobj.addcoupling('mat', 'J2', 'bond', 1);
            % Use rng seed for reproducible results
            swobj.optmagk('kbase', kbase_opts, 'seed', 1);

            expected_k = acos(-J2/(2*J1))/(2*pi);
            rel_tol = 1e-5;
            if abs(expected_k - swobj.mag_str.k(1)) > rel_tol*expected_k
                % If this k doesn't match, try 1-k
                expected_k = 1 - expected_k; % k and 1-k are degenerate
            end
            expected_mag_str = testCase.default_mag_str;
            expected_mag_str.k = [expected_k; expected_k; 0];
            testCase.verify_val(swobj.mag_str, expected_mag_str, ...
                                'rel_tol', 1e-5);
        end
        function test_afm_chain_ndbase_pso_varargin_passed(testCase)
            testCase.swobj.addmatrix('label', 'J1', 'value', 1);
            testCase.swobj.addcoupling('mat', 'J1', 'bond', 1);
            % Verify in default case there is no warning
            testCase.verifyWarningFree(...
                @() testCase.swobj.optmagk, ...
                'pso:convergence');
            % Test that MaxIter gets passed through to ndbase.pso, triggers
            % convergence warning
            testCase.verifyWarning(...
                @() testCase.swobj.optmagk('MaxIter', 1), ...
                'pso:convergence');
        end
    end
end
