classdef unittest_spinw_optmagstr < sw_tests.unit_tests.unittest_super

    properties
        tri = [];
        afc = [];
        opt_tri_mag_str = struct('nExt', int32([1 1 1]), ...
                                 'k', [1/3; 1/3; 0], ...
                                 'F', [1; 1i; 0]);
        tri_optmagstr_args = {'func', @gm_planar, ...
                              'xmin', [0 0 0 0 0 0], ...
                              'xmax', [0 1/2 1/2 0 0 0]};
        orig_rng_state = []
    end
    properties (TestParameter)
        xparams = {'xmin', 'xmax', 'x0'};
        optparams = {{'maxfunevals', 5}, ...
                     {'maxiter', 10}, ...
                     {'tolx', 1e-3}, ...
                     {'tolfun', 1e-4}, ...
                     {'maxfunevals', 5, 'maxiter', 10}};
    end
    methods (Static)
        function [S, k, n] = optmagstr_custom_func(S0, x)
            S = [1; 0; 0];
            k = [1/3 1/3 0];
            n = [0 0 1];
        end
    end
    methods (TestClassSetup)
        function set_seed(testCase)
            testCase.orig_rng_state = rng;
            rng('default');
        end
    end
    methods (TestMethodSetup)
        function setup_afm_tri(testCase)
            testCase.tri = spinw();
            testCase.tri.genlattice('lat_const',[3 3 9],'angled',[90 90 120]);
            testCase.tri.addatom('r',[0 0 0],'S',1);
            testCase.tri.gencoupling('maxDistance',10);
            testCase.tri.addmatrix('value', 1,'label','J1');
            testCase.tri.addcoupling('mat', 'J1','bond', 1);
        end
        function setup_afc(testCase)
            % From tutorial 22
            testCase.afc = spinw();
            testCase.afc.genlattice('lat_const',[3 4 4],'angled',[90 90 90]);
            testCase.afc.addatom('r',[0 0 0],'S',1);
            testCase.afc.addmatrix('label', 'A', 'value', diag([0 0 0.1]));
            testCase.afc.addmatrix('label','J1', 'value', 1);
            testCase.afc.addmatrix('label','J2', 'value', 1/3);
            testCase.afc.gencoupling;
            testCase.afc.addcoupling('mat', 'J1', 'bond', 1);
            testCase.afc.addcoupling('mat', 'J2', 'bond', 5);
            testCase.afc.addaniso('A');
        end
    end
    methods (TestMethodTeardown)
        function reset_seed(testCase)
            rng(testCase.orig_rng_state);
        end
    end
    methods (Test)
        function test_no_mag_atom_throws_error(testCase)
            swobj = spinw();
            testCase.verifyError(@() swobj.optmagstr, 'spinw:optmagstr:NoMagAtom');
        end

        function test_wrong_xparam_length_warns(testCase, xparams)
            params = struct('func', @gm_planar, ...
                            'xmin', [0 0 0 0 0 0], ...
                            'xmax', [0 1/2 1/2 0 0 0], ...
                            'x0', [0 1/4 1/4 0 0 0]);
            xparam = params.(xparams);
            params.(xparams) = xparam(1:end-1);
            testCase.verifyWarning(...
                @() testCase.tri.optmagstr(params), ...
                'spinw:optmagstr:WrongLengthXParam');
        end

        function test_optmagstr_tri_af_out_planar_xmin_xmax(testCase)
            out = testCase.tri.optmagstr(testCase.tri_optmagstr_args{:});
            xmin = testCase.tri_optmagstr_args{4};
            xmax = testCase.tri_optmagstr_args{6};
 
            % Note double {} in 'xname', 'boundary' below, otherwise MATLAB
            % creates a cell array of structs
            expected_out = struct( ...
                'x', [0 1/3 1/3 0 0 0], ...
                'e', -1.5, ...
                'exitflag', 1, ...
                'param', struct('epsilon', 1e-5, ...
                                'func', @gm_planar, ...
                                'boundary', {{'per', 'per', 'per'}}, ...
                                'xmin', xmin, ...
                                'xmax', xmax, ...
                                'x0', [], ...
                                'tolx', 1e-4, ...
                                'tolfun', 1e-5, ...
                                'maxfunevals', 1e7, ...
                                'nRun', 1, ...
                                'maxiter', 1e4, ...
                                'title', 'Optimised magnetic structure using simplex search', ...
                                'tid', 1), ...
                'fname', '2D planar structure', ...
                'xname', {{'Phi1_rad'  'kx_rlu'  'ky_rlu'  'kz_rlu'  'nTheta_rad'  'nPhi_rad'}}, ...
                'title', 'Optimised magnetic structure using simplex search');
            % Some values will change on each call, just check they
            % exist and are of the right type
            assert(isa(out.datestart, 'char'));
            assert(isa(out.dateend, 'char'));
            assert(isa(out.output.iterations, 'double'));
            assert(isa(out.output.funcCount, 'double'));
            assert(isa(out.output.algorithm, 'char'));
            assert(isa(out.output.message, 'char'));
            testCase.verify_obj(out.obj, testCase.tri);
            testCase.verify_val( ...
                rmfield(out, {'output', 'datestart', 'dateend', 'obj'}), ...
                expected_out, 'rel_tol', 1e-3);

            testCase.verify_val(testCase.tri.mag_str, testCase.opt_tri_mag_str, ...
                                'rel_tol', 1e-3);
        end

        function test_optmagstr_tri_af_nExt_init(testCase)
            % Test that if a magnetic structure is initialised with nExt,
            % it is used in optmagstr
            testCase.disable_warnings('spinw:genmagstr:SnParallel');
            testCase.tri.genmagstr('mode', 'random', 'nExt', [3 1 1]);
            testCase.tri.optmagstr('func', @gm_planar, ...
                                   'xmin', [0 pi/2 pi 0 0 0 0 0], ...
                                   'xmax', [0 pi 3*pi/2 0 1/2 0 0 0]);
            testCase.verify_val(testCase.tri.mag_str.k, [0; 1/3; 0], ...
                                   'rel_tol', 1e-3);
        end

        function test_optmagstr_tri_af_named_xparam(testCase)
            % Test using named params
            testCase.tri.optmagstr('func', @gm_planar, ...
                                   'Phi1_rad', [0 0], ...
                                   'kx_rlu', [0 0.5], ...
                                   'ky_rlu', [0 0.5], ...
                                   'kz_rlu', [0 0], ...
                                   'nTheta_rad', [0 0], ...
                                   'nPhi_rad', [0 0]);
            testCase.verify_val(testCase.tri.mag_str, testCase.opt_tri_mag_str, ...
                                'rel_tol', 1e-3);
        end

        function test_optmagstr_tri_af_x0(testCase)
            % Test initialising near a min converges to min
            diff = 0.05;
            testCase.tri.optmagstr('func', @gm_planar, 'x0', [0 1/3+diff 1/3+diff 0 0 0]);
            expected_mag_str =  testCase.opt_tri_mag_str;
            expected_mag_str.F = [-1i; 1; 0];
            % Use abs tol for F = 0
            testCase.verify_val(testCase.tri.mag_str, expected_mag_str, ...
                                'rel_tol', 1e-3, 'abs_tol', 1e-3);
        end

        function test_optmagstr_tri_af_custom_func(testCase)
            testCase.tri.optmagstr('func', @testCase.optmagstr_custom_func, 'xmin', [0], 'xmax', [0]);
            testCase.verify_val(testCase.tri.mag_str, testCase.opt_tri_mag_str, ...
                                'rel_tol', 1e-3);
        end

        function test_optmagstr_tri_af_custom_func_requires_xmin_and_xmax(testCase)
            custom_func = @testCase.optmagstr_custom_func;
            testCase.verifyError(@() testCase.tri.optmagstr('func', custom_func), ...
                                 'spinw:optmagtr:WrongInput');
            testCase.verifyError(@() testCase.tri.optmagstr('func', custom_func, 'xmin', [0]), ...
                                 'spinw:optmagtr:WrongInput');
            testCase.verifyError(@() testCase.tri.optmagstr('func', custom_func, 'xmax', [0]), ...
                                 'spinw:optmagtr:WrongInput');
        end

        function test_optmagstr_tri_af_custom_func_wrong_number_of_outputs(testCase)
            function [S, k] = custom_func(S0, x)
                S = [1; 0; 0];
                k = [0 0 0];
            end
            testCase.verifyError(@() testCase.tri.optmagstr(...
                'func', @custom_func, 'xmin', [0], 'xmax', [0]), ...
                'MATLAB:TooManyOutputs');
        end

        function test_optmagstr_tri_af_custom_func_wrong_number_of_inputs(testCase)
            function [S, k, n] = custom_func(S0)
                S = [1; 0; 0];
                k = [0 0 0];
                n = [0 0 1];
            end
            testCase.verifyError(@() testCase.tri.optmagstr(...
                'func', @custom_func, 'xmin', [0], 'xmax', [0]), ...
                'MATLAB:TooManyInputs');
        end

        function test_optmagstr_tri_af_epsilon(testCase)
            testCase.disable_warnings('spinw:genmagstr:SnParallel');
            % Test that large epsilon doesn't rotate spins
            testCase.tri.optmagstr('epsilon', 1.);
            expected_mag_str = testCase.opt_tri_mag_str;
            expected_mag_str.k = [0 0 0]';
            expected_mag_str.F = [0 0 1]';
            testCase.verify_val(testCase.tri.mag_str, expected_mag_str, ...
                                'rel_tol', 1e-3);
        end

        function test_optmagstr_tri_nRun(testCase)
            % sw_timeit called in each nRun loop, and before and after
            nRun = 4;
            mock_sw_timeit = sw_tests.utilities.mock_function('sw_timeit');
            testCase.tri.optmagstr(testCase.tri_optmagstr_args{:}, 'nRun', nRun);
            testCase.assertEqual(mock_sw_timeit.n_calls, nRun + 2);
            testCase.verify_val(testCase.tri.mag_str, testCase.opt_tri_mag_str, ...
                                'rel_tol', 1e-3);
        end

        function test_optmagstr_tri_set_title(testCase)
            title = 'Test';
            out = testCase.tri.optmagstr(testCase.tri_optmagstr_args{:}, 'title', title);
            testCase.assertEqual(out.title, title);
            testCase.verify_val(testCase.tri.mag_str, testCase.opt_tri_mag_str, ...
                                'rel_tol', 1e-3);
        end

        function test_optmagstr_tri_tid(testCase)
            sw_timeit_mock = sw_tests.utilities.mock_function('sw_timeit');
            tid = 2;
            testCase.tri.optmagstr(testCase.tri_optmagstr_args{:}, 'tid', tid);
            % check tid used in timing
            for irow = 1:sw_timeit_mock.n_calls
                testCase.assertEqual(sw_timeit_mock.arguments{irow}{3}, tid)
            end
            testCase.verify_val(testCase.tri.mag_str, testCase.opt_tri_mag_str, ...
                                'rel_tol', 1e-3);
        end

        function test_optmagstr_optimisation_params(testCase, optparams)
            xmin = [0 0  0 0 0 0 0];
            xmax = [pi/2 0 1/2 0 0 0 0];
            mock_optimset = sw_tests.utilities.mock_function('optimset', ...
                optimset('Display', 'off', optparams{:}));
            testCase.disable_warnings('spinw:genmagstr:SnParallel');
            testCase.afc.optmagstr('xmin', xmin, 'xmax', xmax, ...
                                   optparams{:});
            testCase.assertEqual(mock_optimset.n_calls, 1);
            argslower = cellfun(@(c) lower(c), mock_optimset.arguments{1}(1:2:end), 'UniformOutput', false);
            for ii = find(ismember(argslower, optparams(1:2:end)))
                jj = find(ismember(optparams(1:2:end), argslower{ii}));
                testCase.verifyEqual(mock_optimset.arguments{1}{2*ii}, optparams{2*jj});
            end
        end

        function test_afc_no_init(testCase)
            % Test without initialising parameters, doesn't converge k
            converged_k = [0.385; 0; 0];
            testCase.afc.optmagstr();
            actual_k = testCase.afc.mag_str.k;
            testCase.verifyGreaterThan(sum(abs(converged_k - actual_k)), 0.1);
        end

        function test_afc_gm_spherical3d(testCase)
            xmin = [0 0  0 0 0 0 0];
            xmax = [pi/2 0 1/2 0 0 0 0];
            testCase.afc.optmagstr(...
                'func', @gm_spherical3d, 'xmin', xmin, 'xmax', xmax);
            expected_k = [0.385; 0; 0];
            testCase.verify_val(testCase.afc.mag_str.k, expected_k, ...
                                'rel_tol', 1e-3, 'abs_tol', 1e-3);
        end

        function test_dm_multiatom_spherical3d(testCase)
            sq = spinw();
            sq.genlattice('lat_const', [4 4 4], 'angled', [90 90 90]);
            sq.addatom('r', [0 0 0], 'S', 1);
            sq.addatom('r', [0.5 0.5 0.5], 'S', 1);
            sq.gencoupling('maxDistance', 10);
            % This is the DM interaction, with the vector along [1 1 1]
            sq.addmatrix('value', [1 1 1], 'label', 'DM');
            sq.addcoupling('mat', 'DM', 'bond', 1);
            sq.addmatrix('value', 1, 'label', 'J1');
            sq.addcoupling('mat', 'J1', 'bond', 2);
            testCase.disable_warnings('spinw:genmagstr:SnParallel');
            % Sometimes fails to find min, run multiple times
            sq.optmagstr('func', @gm_spherical3d, ...
                         'xmin', [-pi/2 -pi -pi/2 -pi, 0 0 0, 0 0], ...
                         'xmax', [pi/2 pi pi/2 pi, 0 0 0, 0 0], ...
                         'nRun', 5);
            spin_angles = {{90, 121}, {31, 180}}; % theta, phi
            expected_F = zeros(3, length(spin_angles));
            for i=1:length(spin_angles)
                [theta, phi] = spin_angles{i}{:};
                expected_F(1, i) = sind(theta)*cosd(phi); % a
                expected_F(2, i) = sind(theta)*sind(phi); % b
                expected_F(3, i) = cosd(theta); % c
            end
            testCase.verify_val(sq.mag_str.F, expected_F, 'abs_tol', 0.02);
        end

    end

end
