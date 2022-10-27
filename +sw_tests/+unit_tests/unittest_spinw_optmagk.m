classdef unittest_spinw_optmagk < sw_tests.unit_tests.unittest_super
    % Runs through unit test for @spinw/spinwave.m

    properties
        swobj = [];
        swobj_tri = [];
        % output from optmagk
        default_mag_str = struct('nExt', int32([1 1 1]), ...
                                 'k', [1; 0; 0]);
    end
    methods (TestClassSetup)
        function setup_chain_model(testCase)            
            testCase.swobj = spinw();
            testCase.swobj.genlattice('lat_const', [3 8 8])
            testCase.swobj.addatom('r',[0; 0; 0],'S',1)
            testCase.swobj.gencoupling();
        end
    end
    methods (Test)
        function test_fm_opt(testCase)
            testCase.swobj.addmatrix('label', 'J1', 'value', -1);
            testCase.swobj.addcoupling('mat', 'J1', 'bond', 1);
            out = testCase.swobj.optmagk;
            out.stat = rmfield(out.stat, 'nFunEvals');
            k = [1; 0; 0];
            F = [sqrt(1/3) + 1i*sqrt(1/2); ...
                 sqrt(1/3); ...
                 sqrt(1/3) - 1i*sqrt(1/2)];
            expected_out = struct('k', k, ...
                                  'E', -1, ....
                                  'F', F, ...
                                  'stat', struct('S', 0, ...
                                                 'exitflag', -1));
            % Test struct output by optmagk
            testCase.verify_val(out, expected_out, 'abs_tol', 1e-4);
            % Also test spinw attributes have been set
            expected_mag_str = testCase.default_mag_str;
            expected_mag_str.F = F;
            testCase.verify_val(testCase.swobj.mag_str, expected_mag_str, ...
                                'abs_tol', 1e-4);
        end
    end
    methods (Test)
        function test_afm_opt(testCase)
            testCase.swobj.addmatrix('label', 'J1', 'value', 1);
            testCase.swobj.addcoupling('mat', 'J1', 'bond', 1);
            testCase.swobj.optmagk;
            % Also test spinw attributes have been set
            expected_mag_str = testCase.default_mag_str;
            expected_mag_str.k = [0.5; 0; 0];
            expected_mag_str.F = [sqrt(1/3) + 1i*sqrt(1/2); ...
                                  sqrt(1/3); ...
                                  sqrt(1/3) - 1i*sqrt(1/2)];
            testCase.verify_val(testCase.swobj.mag_str, expected_mag_str, ...
                                'abs_tol', 1e-4);
        end
    end
end
