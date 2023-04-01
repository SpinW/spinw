classdef unittest_spinw_optmagsteep < sw_tests.unit_tests.unittest_super

    properties
        swobj = [];
        default_magstr = struct('S', [0  0; 0 0; -1 1], 'n', [0 0 1], ...
                                'N_ext', [2 1 1], 'k', [0 0 0], ...
                                'exact', true)
        % output from optmagsteep
        default_opt = struct('M', [0 0; 0 0; -1 1],...
                             'dM', 6.6e-17,...
                             'e', -1.1, ...
                             'nRun', 1, ...
                             'title', ['Optimised magnetic structure ', ...
                                      'using the method of steepest ', ...
                                      'descent']);
       orig_rng_state = []
    end
    properties (TestParameter)
        existing_plot = {true, false};
    end
    methods (TestClassSetup)
        function set_seed(testCase)
            testCase.orig_rng_state = rng;
            rng('default');
        end
    end
    methods (TestMethodSetup)
        function setup_afm_chain_model_easy_axis(testCase)            
            testCase.swobj = spinw();
            testCase.swobj.genlattice('lat_const', [2 3 6])
            testCase.swobj.addatom('r',[0; 0; 0],'S',1)
            testCase.swobj.addmatrix('label', 'A', ...
                                     'value', diag([0 0 -0.1])) % c easy
            testCase.swobj.addmatrix('label', 'J1', 'value', 1)  % AFM
            testCase.swobj.gencoupling('maxDistance', 6);
            testCase.swobj.addcoupling('mat', 'J1', 'bond', 1); % along a
            testCase.swobj.addaniso('A');
        end
    end
    methods(TestMethodTeardown)
        function reset_seed(testCase)
            rng(testCase.orig_rng_state);
        end
    end
    methods (Test)
        function test_only_one_spin_throws_error(testCase)
            testCase.swobj.genmagstr('mode', 'direct', 'S', [0; 0; 1], ...
                                     'k',[0, 0, 0]);
            func_call = @() testCase.swobj.optmagsteep('random', true);
            testCase.verifyError(func_call, 'spinw:optmagsteep:NoField')
        end
        
        function test_invalid_nExt_throws_error(testCase)
            func_call = @() testCase.swobj.optmagsteep('nExt', [0, 0, 0]);
            testCase.verifyError(func_call, 'spinw:optmagsteep:WrongInput')
        end
        
        function test_warns_if_not_converged(testCase)
            % init moment along hard-axis (far from minima) and run 1 iter
            testCase.swobj.genmagstr('mode', 'helical', 'S', [1; 0; 0], ...
                                     'k',[0.5,0,0], 'n', [0,1,0], ...
                                     'nExt', [2,1,1]);
            func_call = @() testCase.swobj.optmagsteep('nRun', 1);
            testCase.verifyWarning(func_call, ...
                                   'spinw:optmagsteep:NotConverged')
        end
        
        function test_warns_of_self_coupling_spins(testCase)
            % add AFM exchange to bond along b (will double cell along b)
            testCase.swobj.addmatrix('label', 'J2', 'value', 0.25)  % AFM
            testCase.swobj.addcoupling('mat', 'J2', 'bond', 2); % along b
            % run opt with k only along a
            func_call = @() testCase.swobj.optmagsteep('NExt', [2, 1, 1], ...
                                                       'nRun', 1);
            testCase.verifyWarning(func_call, ...
                                   'spinw:optmagsteep:SelfCoupling')
        end
        
        function test_starting_from_groundstate_executes_one_iteration(testCase)
            % generate ground state magnetic structure
            testCase.swobj.genmagstr('mode', 'helical', 'S', [0; 0; -1], ...
                                     'k',[0.5,0,0], 'n', [0,1,0], ...
                                     'nExt', [2,1,1]);
            opt = testCase.swobj.optmagsteep;
            testCase.verify_val(testCase.swobj.magstr, ...
                                testCase.default_magstr, 'abs_tol', 1e-6);
            fields = {'datestart', 'dateend', 'obj', 'param'};
            testCase.verify_val(rmfield(opt, fields), testCase.default_opt);
            testCase.verify_obj(opt.obj, testCase.swobj)
        end
        
        function test_converges_local_minimum_along_hard_axis(testCase)
            % FM with S//a (hard axis) - no component along easy-axis
            testCase.swobj.genmagstr('mode', 'direct', ...
                                     'S', [1 1; 0 0; 0 0], ...
                                     'k',[0,0,0], 'nExt', [2,1,1]);
            testCase.swobj.optmagsteep;  % results in AFM S//a
            expected_magstr = testCase.default_magstr;
            expected_magstr.S = [-1, 1; 0, 0; 0, 0];
            testCase.verify_val(testCase.swobj.magstr, expected_magstr);
            testCase.verify_val(testCase.swobj.energy, -1.0)
        end
        
        function test_spins_align_along_easy_axis_external_field(testCase)
            testCase.swobj.field([0 0 1]);
            % FM with S//a (hard axis) - no component along easy-axis
            testCase.swobj.genmagstr('mode', 'direct', ...
                                     'S', [1 1; 0 0; 0 0], ...
                                     'k',[0,0,0], 'nExt', [2,1,1]);
            testCase.swobj.optmagsteep('nRun', 150);  % results in AFM S//c
            expected_magstr = testCase.default_magstr;
            expected_magstr.S = [0, 0; 0, 0; 1, -1];
            testCase.verify_val(testCase.swobj.magstr, expected_magstr, ...
                'abs_tol', 1e-6);
            testCase.verify_val(testCase.swobj.energy, -1.1, ...
                'abs_tol', 1e-10)  % same energy as grd state without field
        end
        
        function test_converges_global_minimum(testCase)
            % FM with component of S // c (easy-axis)
            testCase.swobj.genmagstr('mode', 'direct', ...
                                     'S', [0 0; 1 1; 1 1], ...
                                     'k',[0,0,0], 'nExt', [2,1,1]);
            % try with one iterations - will not convergence
            testCase.verifyWarning(...
                @() testCase.swobj.optmagsteep('nRun', 1), ...
                'spinw:optmagsteep:NotConverged');
            testCase.swobj.optmagsteep('nRun', 150);  % results in AFM S//c
            testCase.verify_val(testCase.swobj.magstr, ...
                                testCase.default_magstr, 'abs_tol', 1e-8);
            testCase.verify_val(testCase.swobj.energy, -1.1)
        end
        
        function test_random_init_of_spins(testCase)
            % FM with S//a (hard axis) - no component along easy-axis
            % optmagsteep would converge at local not global minimum
            testCase.swobj.genmagstr('mode', 'direct', ...
                                     'S', [1 1; 0 0; 0 0], ...
                                     'k',[0,0,0], 'nExt', [2,1,1]);
            testCase.swobj.optmagsteep('random', true, 'nRun', 200)
            expected_magstr = testCase.default_magstr;
            % first moment can be up/down (same energy) so adjust 
            % expected value to have same z-component sign on 1st S
            expected_magstr.S = -sign(testCase.swobj.magstr.S(end,1))*...
                expected_magstr.S;
            testCase.verify_val(testCase.swobj.magstr, ...
                                expected_magstr, 'abs_tol', 1e-6);
            testCase.verify_val(testCase.swobj.energy, -1.1);
        end
        
        function test_random_init_spins_if_no_initial_magstr(testCase)
            testCase.swobj.optmagsteep('nExt', [2,1,1], 'nRun', 250);
            testCase.verify_val(testCase.swobj.energy, -1.1);
        end
        
        function test_output_to_fid(testCase)
           % generate ground state magnetic structure
            testCase.swobj.genmagstr('mode', 'helical', 'S', [0; 0; 1], ...
                                     'k',[0.5,0,0], 'n', [0,1,0], ...
                                     'nExt', [2,1,1]);
            mock_fprintf = sw_tests.utilities.mock_function('fprintf0');
            fid = 3;
            testCase.swobj.optmagsteep('fid', fid)
            testCase.assertEqual(mock_fprintf.n_calls, 2);
            % check fid used to write file
            for irow = 1:mock_fprintf.n_calls
                testCase.assertEqual(mock_fprintf.arguments{irow}{1}, fid)
            end
        end
        
        function test_plot_moment_each_iteration(testCase, existing_plot)
            testCase.disable_warnings('MATLAB:dispatcher:nameConflict');
            if existing_plot
                existing_fig = testCase.swobj.plot();
            end
            % generate ground state magnetic structure
            testCase.swobj.genmagstr('mode', 'helical', 'S', [0; 0; -1], ...
                                     'k',[0.5,0,0], 'n', [0,1,0], ...
                                     'nExt', [2,1,1]);
            mock_pause = sw_tests.utilities.mock_function('pause');
            tpause = 1e-5;
            testCase.swobj.optmagsteep('plot', true, 'pause', tpause);
            % check magstr plotted (with arrow)
            fig = swplot.activefigure;
            testCase.assertEqual(numel(findobj(fig,'Tag', 'arrow')), 2)
            if existing_plot
                testCase.assertEqual(fig, existing_fig);
            end
            close(fig);
            testCase.assertEqual(mock_pause.n_calls, 1);
            testCase.assertEqual(mock_pause.arguments{1}, {tpause})
        end
        
        function test_convergence_stops_with_given_xtol(testCase)
            testCase.swobj.genmagstr('mode', 'direct', ...
                         'S', [0 0; 1 1; 1 1], ...
                         'k',[0,0,0], 'nExt', [2,1,1]);
            dM_tol = 1e-3;
            opt_struct = testCase.swobj.optmagsteep('TolX', dM_tol, ...
                                                    'saveAll', true);
            testCase.verify_val(opt_struct.dM, dM_tol, 'abs_tol', 1e-4);
            % check M saved for each iteration
            testCase.assertEqual(size(opt_struct.M, 3), opt_struct.nRun);
        end
        
        function test_not_move_moments_in_field_more_than_Hmin(testCase)
            testCase.swobj.genmagstr('mode', 'direct', ...
                         'S', [0 0; 1 1; 1 1], ...
                         'k',[0,0,0], 'nExt', [2,1,1]);
            opt_struct = testCase.swobj.optmagsteep('Hmin', 2);
            % check moments haven't moved
            testCase.verify_val(opt_struct.dM, 0);
            % check structure is unchanged
            expected_magstr = testCase.default_magstr;
            expected_magstr.S = [0, 0; 1, 1; 1, 1]/sqrt(2);
            testCase.verify_val(testCase.swobj.magstr, expected_magstr);
        end
        
        function test_multiple_atoms_in_unit_cell(testCase)
            testCase.swobj.addatom('r',[0.5; 0.5; 0.5],'S',1)
            testCase.swobj.gencoupling('maxDistance', 6, 'dMin', 0.1);
            testCase.swobj.addaniso('A'); % add again as cleared above
            % add AFM coupling of spins in same unit cell
            testCase.swobj.addcoupling('mat', 'J1', 'bond', 3);
            testCase.swobj.genmagstr('mode', 'direct', 'S', [0 0; 0 0; 1 1], ...
                                    'k',[0, 0, 0]); % FM initial state
                                
            testCase.swobj.optmagsteep();
            
            expected_magstr = testCase.default_magstr;
            expected_magstr.N_ext = [1 1 1];
            testCase.verify_val(testCase.swobj.magstr, expected_magstr);
            testCase.verify_val(testCase.swobj.energy, -4.1)
        end
        
    end

end