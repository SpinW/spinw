classdef unittest_spinw_optmagsteep < sw_tests.unit_tests.unittest_super

    properties
        swobj = [];
        default_mag_str = struct('nExt', int32([2 2 1]), ...
                                 'k', [0; 0; 0], ...
                                 'F', repmat([0,        0; 
                                              0,        0; 
                                              -1+0j, 1+0j], 1, 2));
    end
    properties (TestParameter)

    end
    methods (TestClassSetup)
        function setup_afm_chain_model_easy_axis(testCase)            
            testCase.swobj = spinw();
            testCase.swobj.genlattice('lat_const',[2 4 6])
            testCase.swobj.addatom('r',[0; 0; 0],'S',1)
            testCase.swobj.addmatrix('label', 'A', ...
                                     'value', diag([0 0 -0.1])) % c easy
            testCase.swobj.addmatrix('label', 'J1', 'value', 1)  % AFM
            testCase.swobj.gencoupling();
            testCase.swobj.addcoupling('mat', 'J1', 'bond', 1); % along a
            testCase.swobj.addaniso('A');

        end
    end

    methods (Test)
        function test_invalid_nExt_throws_error(testCase)
            testCase.swobj.genmagstr('mode', 'direct', 'S', [0; 0; 1], ...
                                     'k',[0,0,0]);
            func_call = @() testCase.swobj.optmagsteep('random', true);
            testCase.verifyError(func_call, 'spinw:optmagsteep:NoField')
        end
        
        function test_converges_local_minimum_along_hard_axis(testCase)
            % FM with S//a (hard axis) - no component along easy-axis
            testCase.swobj.genmagstr('mode', 'direct', ...
                                     'S', repmat([1;0;0], 1,4), ...
                                     'k',[0,0,0], 'nExt', [2,2,1]);
            testCase.swobj.optmagsteep;  % results in AFM S//a
            expected_magstr = testCase.default_mag_str;
            expected_magstr.F = repmat([-1, 1; -1j, 1j; 0, 0], 1, 2);
            testCase.verify_val(testCase.swobj.mag_str, expected_magstr);
            testCase.verify_val(testCase.swobj.energy, -1.0)
        end
        
        function test_converges_global_minimum(testCase)
            % FM with component of S // c (easy-axis)
            testCase.swobj.genmagstr('mode', 'direct', ...
                                     'S', repmat([0;1;1], 1,4), ...
                                     'k',[0,0,0], 'nExt', [2,2,1]);
            % try with one iterations - will not convergence
            testCase.verifyWarning(...
                @() testCase.swobj.optmagsteep('nRun', 1), ...
                'spinw:optmagsteep:NotConverged');
            testCase.swobj.optmagsteep('nRun', 200);  % results in AFM S//c
            expected_magstr = testCase.default_mag_str;
            expected_magstr.F = repmat([0, 0; 0, 0; -1+0j, 1+0j], 1, 2);
            testCase.verify_val(testCase.swobj.mag_str, expected_magstr, ...
                                'abs_tol', 1e-8);
            
            testCase.verify_val(testCase.swobj.energy, -1.1)
        end
        
        function test_random_init_of_spins(testCase)
            % FM with S//a (hard axis) - no component along easy-axis
            testCase.swobj.genmagstr('mode', 'direct', ...
                                     'S', repmat([1;0;0], 1,4), ...
                                     'k',[0,0,0], 'nExt', [2,2,1]);
            has_easy_comp = false;
            for irun = 1:5
                testCase.swobj.optmagsteep('random', true)
                has_easy_comp = abs(real(testCase.swobj.mag_str.F(end,1))) > 0;
                if has_easy_comp
                    break;
                end
            end
            if has_easy_comp
                expected_magstr = testCase.default_mag_str;
                expected_magstr.F = repmat([0, 0; 0, 0; 1+0j, 1+0j], ...
                                           1, 2);
                % no coupling between chains so can be AFM and FM ordered
                % along b-axis
                expected_magstr.F(end,:) = expected_magstr.F(end,:).*sign(...
                    testCase.swobj.mag_str.F(end,:));
                testCase.verify_val(testCase.swobj.mag_str, ...
                                    expected_magstr, 'abs_tol', 5e-8);
                testCase.verify_val(testCase.swobj.energy, -1.1)
            end
        end
        
    end

end