classdef unittest_spinw_optmagk < sw_tests.unit_tests.unittest_super
    % Runs through unit test for @spinw/spinwave.m

    properties
        swobj = [];
        % output from optmagk
        default_mag_str = struct('nExt', int32([1 1 1]), ...
                                 'F', [sqrt(1/3) + 1i*sqrt(1/2); ...
                                       sqrt(1/3); ...
                                       sqrt(1/3) - 1i*sqrt(1/2)], ...
                                 'k', [1; 0; 0]);
    end
    properties (TestParameter)
        kbase_opts = {[1; 1; 0], [1 0; 0 1; 0 0]};
    end
    methods (TestMethodSetup)
        function setup_chain_model(testCase)            
            testCase.swobj = spinw();
            testCase.swobj.genlattice('lat_const', [3 8 8])
            testCase.swobj.addatom('r',[0; 0; 0],'S',1)
            testCase.swobj.gencoupling();
        end
    end
    methods (Test)
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
        function test_fm_chain_optk(testCase)
            testCase.swobj.addmatrix('label', 'J1', 'value', -1);
            testCase.swobj.addcoupling('mat', 'J1', 'bond', 1);
            out = testCase.swobj.optmagk;
            out.stat = rmfield(out.stat, 'nFunEvals');

            expected_mag_str = testCase.default_mag_str;
            expected_out = struct('k', expected_mag_str.k, ...
                                  'E', -1, ....
                                  'F', expected_mag_str.F, ...
                                  'stat', struct('S', 0, ...
                                                 'exitflag', -1));
            % Test struct output by optmagk
            testCase.verify_val(out, expected_out, 'abs_tol', 1e-4);
            % Also test spinw attributes have been set
            expected_mag_str = testCase.default_mag_str;
            testCase.verify_val(testCase.swobj.mag_str, expected_mag_str, ...
                                'abs_tol', 1e-4);
        end
        function test_afm_chain_optk(testCase)
            testCase.swobj.addmatrix('label', 'J1', 'value', 1);
            testCase.swobj.addcoupling('mat', 'J1', 'bond', 1);
            testCase.swobj.optmagk;
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
            swobj.optmagk('kbase', kbase_opts);

            expected_k = acos(-J2/(2*J1))/(2*pi);
            rel_tol = 1e-5;
            if abs(expected_k - swobj.mag_str.k(1)) > rel_tol*expected_k
                % If this k doesn't match, try 1-k
                expected_k = 1 - expected_k; % k and 1-k are degenerate
            end
            expected_mag_str = testCase.default_mag_str;
            expected_mag_str.k = [expected_k; expected_k; 0];
            testCase.verify_val(swobj.mag_str, expected_mag_str, ...
                                'rel_tol', rel_tol);
        end
    end
end
