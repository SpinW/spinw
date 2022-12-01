classdef unittest_spinw_optmagstr < sw_tests.unit_tests.unittest_super

    properties
        tri = [];
        opt_tri_mag_str = struct('nExt', int32([1 1 1]), ...
                                 'k', [1/3; 1/3; 0], ...
                                 'F', [1; 1i; 0]);
        tri_optmagstr_args = {'func', @gm_planar, ...
                              'xmin', [0 0 0 0 0 0], ...
                              'xmax', [0 1/2 1/2 0 0 0]};
       orig_rng_state = []
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
            testCase.tri.genmagstr('mode','helical','S',[1 0 0]','k',[0 0 0]);
        end
    end
    methods(TestMethodTeardown)
        function reset_seed(testCase)
            rng(testCase.orig_rng_state);
        end
    end
    methods (Test)
        function test_no_mag_atom_throws_error(testCase)
            swobj = spinw();
            testCase.verifyError(@() swobj.optmagstr, 'spinw:optmagstr:NoMagAtom');
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

        function test_optmagstr_tri_af_no_init(testCase)
            % Test without initialising parameters, converges to different
            % k
            testCase.tri.optmagstr();
            converged_k = [1/3; 1/3; 0];
            actual_k = testCase.tri.mag_str.k;
            testCase.verifyGreaterThan(sum(abs(converged_k - actual_k)), 0.1);
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

        function test_optmagstr_tri_af_epsilon(testCase)
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

        function test_optmagstr_optimisation_params(testCase)
            % We could just mock optimset and check correct args are passed
            % through, but Matlab doesn't allow mocking functions so just
            % check some behaviour and the output struct. Our mock_function
            % doesn't support complex return values such as structs.
            tolx = 1e-5;
            tolfun = 1e-6;
            maxfunevals = 3;
            maxiter = 2;
            out = testCase.tri.optmagstr(testCase.tri_optmagstr_args{:}, ...
                                         'tolx', tolx, 'tolfun', tolfun, ...
                                         'maxfunevals', maxfunevals, ...
                                         'maxiter', maxiter);
            converged_k = [1/3; 1/3; 0];
            actual_k = testCase.tri.mag_str.k;
            % Test with low iterations k doesn't converge
            testCase.verifyGreaterThan(sum(abs(converged_k - actual_k)), 0.1);
            % Test that params are in output struct
            testCase.verifyEqual(out.param.tolx, tolx);
            testCase.verifyEqual(out.param.tolfun, tolfun);
            testCase.verifyEqual(out.param.maxfunevals, maxfunevals);
            testCase.verifyEqual(out.param.maxiter, maxiter);
        end

    end

end