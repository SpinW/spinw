classdef unittest_spinw_optmagstr < sw_tests.unit_tests.unittest_super

    properties
        swobj = [];
        tri = [];
        afm = [];
        afm3d = [];
        fr_afm = [];
        default_tri_mag_str = struct('nExt', int32([1 1 1]), ...
                                 'k', [1/3; 1/3; 0], ...
                                 'F', [1; 1i; 0]);
        tri_optmagstr_args = {'func', @gm_planar, ...
                              'xmin', [0 0 0 0 0 0], ...
                              'xmax', [0 1/2 1/2 0 0 0]};
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
            % k = [0 0 1/2] global
            % k = [1/2 0 0] local
            testCase.afm = spinw();
            testCase.afm.genlattice('lat_const', [2 3 6])
            testCase.afm.addatom('r',[0; 0; 0],'S',1)
            testCase.afm.addmatrix('label', 'A', ...
                                   'value', diag([0 0 -0.1])) % c easy
            testCase.afm.addmatrix('label', 'J1', 'value', 1)  % AFM
            testCase.afm.gencoupling('maxDistance', 6);
            testCase.afm.addcoupling('mat', 'J1', 'bond', 1); % along a
            testCase.afm.addaniso('A');
            testCase.afm.genmagstr('mode', 'helical', 'S', [1; 0; 0], ...
                                     'k',[0.5,0,0], 'n', [0,1,0], ...
                                     'nExt', [2,1,1]);
        end
        function setup_frustrated_afm_chain(testCase)
            % Tutorial 3
            % k = [0.2301 0 0]
            testCase.fr_afm = spinw();
            testCase.fr_afm.genlattice('lat_const', [3 8 10])
            testCase.fr_afm.addatom('r',[0; 0; 0], 'S', 1)
            testCase.fr_afm.addmatrix('label', 'J1', 'value', -1)
            testCase.fr_afm.addmatrix('label', 'J2', 'value', 2)
            testCase.fr_afm.gencoupling('maxDistance', 7);
            testCase.fr_afm.addcoupling('mat', 'J1', 'bond', 1);
            testCase.fr_afm.addcoupling('mat', 'J2', 'bond', 2);
            testCase.fr_afm.genmagstr('mode','helical', 'k',[0.25 0 0], ...
                                     'n',[0 0 1], 'S',[1; 0; 0])
        end
        function setup_afm_chain_3d(testCase)
            % Tutorial 22 - maybe dodgy
            testCase.afm3d = spinw();
            testCase.afm3d.genlattice('lat_const', [3 4 4]);
            testCase.afm3d.addatom('r',[0; 0; 0], 'S', 1)
            testCase.afm3d.addmatrix('label', 'A', 'value', diag([0 0 0.1]));
            testCase.afm3d.addmatrix('label', 'J1', 'value', 1);
            testCase.afm3d.addmatrix('label', 'J2', 'value', 1/3);
            testCase.afm3d.gencoupling();
            testCase.afm3d.addcoupling('mat', 'J1', 'bond', 1);
            testCase.afm3d.addcoupling('mat', 'J2', 'bond', 5);
            testCase.afm3d.addaniso('A');
            testCase.afm3d.genmagstr('mode','helical', 'k',[0.25 0 0], ...
                                     'n',[0 0 1], 'S',[1; 0; 0]);
        end
        function setup_afm_tri(testCase)
            testCase.tri = spinw();
            testCase.tri.genlattice('lat_const',[3 3 9],'angled',[90 90 120]);
            testCase.tri.addatom('r',[0 0 0],'S',1);
            testCase.tri.gencoupling('maxDistance',10);
            testCase.tri.addmatrix('value', 1,'label','J1');
            testCase.tri.addcoupling('mat', 1,'bond', 1);
            testCase.tri.genmagstr('mode','helical','S',[1 0 0]','k',[0 0 0]);
            %testCase.tri.optmagstr('func',@gm_planar,'xmin',[0 0 0 0 0 0],'xmax',[0 1/2 1/2 0 0 0],'nRun',10)
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

            testCase.verify_val(testCase.tri.mag_str, testCase.default_tri_mag_str, ...
                                'rel_tol', 1e-3);
        end

        function test_optmagstr_tri_af_no_init(testCase)
            % Test without initialising parameters, a different min is found
            testCase.tri.optmagstr();
            expected_mag_str =  testCase.default_tri_mag_str;
            expected_mag_str.k = [2/3; 2/3; 0.861806743515228];
            expected_mag_str.F = [ 0.675784235861410 + 0.354449639245427i; ...
                                  -0.031366423410885 - 0.862168724764435i; ...
                                  -0.736431812215916 + 0.361981412887755i];
            testCase.verify_val(testCase.tri.mag_str, expected_mag_str, ...
                                'rel_tol', 1e-3);
        end

        function test_optmagstr_tri_af_x0(testCase)
            % Test initialising near a min converges to min
            diff = 0.05;
            testCase.tri.optmagstr('func', @gm_planar, 'x0', [0 1/3+diff 1/3+diff 0 0 0]);
            expected_mag_str =  testCase.default_tri_mag_str;
            expected_mag_str.F = [-1i; 1; 0];
            % Use abs tol for F = 0
            testCase.verify_val(testCase.tri.mag_str, expected_mag_str, ...
                                'rel_tol', 1e-3, 'abs_tol', 1e-3);
        end

        function test_optmagstr_tri_af_epsilon(testCase)
            % Test that large epsilon doesn't rotate spins
            testCase.tri.optmagstr('epsilon', 1.);
            expected_mag_str = testCase.default_tri_mag_str;
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
            testCase.verify_val(testCase.tri.mag_str, testCase.default_tri_mag_str, ...
                                'rel_tol', 1e-3);
        end

        function test_optmagstr_tri_set_title(testCase)
            title = 'Test';
            out = testCase.tri.optmagstr(testCase.tri_optmagstr_args{:}, 'title', title);
            testCase.assertEqual(out.title, title);
            testCase.verify_val(testCase.tri.mag_str, testCase.default_tri_mag_str, ...
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
            testCase.verify_val(testCase.tri.mag_str, testCase.default_tri_mag_str, ...
                                'rel_tol', 1e-3);
        end

end