classdef unittest_spinw_genmagstr < sw_tests.unit_tests.unittest_super

    properties
        swobj = [];
        swobj_tri = [];
        default_mag_str = struct('nExt', int32([1 1 1]), ...
                                 'k', [0; 0; 0], ...
                                 'F', [0; 1; 0]);
    end
    properties (TestParameter)
         fm_chain_input_errors = { ...
             % varargin, identifier
             {{'mode', 'something'}, 'spinw:genmagstr:WrongMode'}; ...
             {{'unit', 'something'}, 'spinw:genmagstr:WrongInput'}; ...
             % nExt can't be zero
             {{'nExt', [0 1 1]}, 'spinw:genmagstr:WrongInput'}; ...
             % n must have dimensions nK, 3
             {{'n', ones(2, 3), 'k', [0 0 1/2]}, 'sw_readparam:ParameterSizeMismatch'}; ...
             % S with direct must have same number of spins as atoms in nExt supercell
             {{'mode', 'direct', 'S', [0; 1; 0], 'nExt', [2 1 1]}, 'spinw:genmagstr:WrongSpinSize'}; ...
             % S with helical must have 1 spin, or same number of spins as
             % atoms in unit or supercell
             {{'mode', 'helical', 'S', [0 1; 1 0; 0 0]}, 'spinw:genmagstr:WrongNumberSpin'}; ...
             {{'mode', 'direct', 'S', [1 0 0]}, 'spinw:genmagstr:WrongInput'}; ...
             {{'mode', 'direct', 'k', [1/2; 0; 0]}, 'spinw:genmagstr:WrongInput'}; ...
             {{'mode', 'direct', 'S', [1 0 0], 'k', [1/2; 0; 0]}, 'spinw:genmagstr:WrongInput'}; ...
             {{'mode', 'tile', 'S', [0 0; 1 1; 1 1]}, 'spinw:genmagstr:WrongInput'}; ...
             % S with helical must be real
             {{'mode', 'helical', 'S', [1.5; 0; 1.5i]}, 'spinw:genmagstr:WrongInput'}; ...
             % Rotate mode must first initialise a magnetic structure
             {{'mode', 'rotate'}, 'spinw:genmagstr:WrongInput'}
             };
         rotate_input_errors = { ...
             % varargin, identifier
             % rotation angle must be real (previously phi=i had special
             % behaviour)
             {{'mode', 'rotate', 'phi', i}, 'spinw:genmagstr:ComplexPhi'}; ...
             % If no angle is supplied to rotate, the rotation axis is set
             % orthogonal to both the first spin and n. If the first spin
             % and n are parallel, this should therefore cause an error
             {{'mode', 'rotate'}, 'spinw:genmagstr:InvalidRotation'}; ...
             };
         ignored_inputs = { ...
             % arguments
             {'mode', 'random', 'S', [1; 0; 0], 'epsilon', 1}, ...
             {'n', [0 1 1], 'mode', 'direct', 'S', [1; 0; 0]}, ...
             {'mode', 'tile', 'k', [0 0 1/2], 'n', [0 0 1], 'S', [1; 0; 0]}, ...
             {'mode', 'helical', 'k', [0 0 1/3], 'S', [1; 0; 0;], 'x0', []}, ...
             {'mode', 'func', 'x0', [pi/2 -pi/4 0 0 [0 0 1/3] pi pi/2], 'unit', 'lu', 'next', [1 2 1]}, ...
             {'mode', 'fourier', 'S', [1; i; 0], 'n', [1 0 0]}
             };
         complex_n_input = {[1 1+i 0], [i 0 0]};
         rotate_phi_inputs = {{'phi', pi/4}, {'phid', 45}, {'phi', pi/4, 'phid', 90}};
         input_norm_output_F = {
             % norm_bool, mag_str.F
             {true, [2; 0; -2i]}, ...
             {false, [1; 0; -1i]}
         };
    end
    methods (TestMethodSetup)
        function setup_chain_model(testCase)
            % Create a simple FM 1D chain model
            testCase.swobj = spinw();
            testCase.swobj.genlattice('lat_const', [3 8 8], 'angled', [90 90 90]);
            testCase.swobj.addatom('r', [0 0 0],'S', 1);
        end
        function setup_tri_model(testCase)
            % Create a simple triangular lattice model
            testCase.swobj_tri = spinw;
            testCase.swobj_tri.genlattice('lat_const', [4 4 6], 'angled', [90 90 120]);
            testCase.swobj_tri.addatom('r', [0 0 0], 'S', 3/2);
        end
    end

    methods (Test)
        function test_invalid_fm_chain_input_raises_error(testCase, fm_chain_input_errors)
            varargin = fm_chain_input_errors{1};
            identifier = fm_chain_input_errors{2};
            testCase.verifyError(...
                @() testCase.swobj.genmagstr(varargin{:}), ...
                identifier)
        end
        function test_invalid_rotate_input_raises_error(testCase, rotate_input_errors)
            varargin = rotate_input_errors{1};
            identifier = rotate_input_errors{2};
            swobj = copy(testCase.swobj);
            swobj.genmagstr('mode', 'direct', 'S', [0; 0; 1]);
            testCase.verifyError(...
                @() swobj.genmagstr(varargin{:}), ...
                identifier)
        end
        function test_no_magnetic_atoms_raises_error(testCase)
            swobj = spinw;
            testCase.verifyError(...
                @() swobj.genmagstr(), ...
                'spinw:genmagstr:NoMagAtom')
        end
        function test_complex_n_input_raises_error(testCase, complex_n_input)
            testCase.verifyError(...
                @() testCase.swobj.genmagstr('mode', 'helical', 'S', [1; 0; 0], 'n', complex_n_input), ...
                'spinw:genmagstr:WrongInput')
        end
        function test_tile_too_few_S_raises_error(testCase)
            swobj = copy(testCase.swobj);
            swobj.addatom('r', [0.5 0.5 0], 'S', 1);
            testCase.verifyError(...
                @() swobj.genmagstr('mode', 'tile', 'S', [0; 1; 1]), ...
                'spinw:genmagstr:WrongInput')
        end
        function test_ignored_input_warns(testCase, ignored_inputs)
            testCase.verifyWarning(...
                @() testCase.swobj.genmagstr(ignored_inputs{:}), ...
                'spinw:genmagstr:UnreadInput')
        end
        function test_rotate_ignored_input_warns(testCase)
            swobj = copy(testCase.swobj);
            % Need to initialise structure before rotating it
            swobj.genmagstr('mode', 'direct', 'S', [1; 0; 0]);
            testCase.verifyWarning(...
                @() swobj.genmagstr('mode', 'rotate', 'k', [0 0 1/3]), ...
                'spinw:genmagstr:UnreadInput')
        end
        function test_helical_spin_size_incomm_with_nExt_warns(testCase)
            swobj = copy(testCase.swobj);
            nExt = [2 1 1];
            k = [1/3 0 0];
            testCase.verifyWarning(...
                @() swobj.genmagstr('mode', 'helical', ...
                                    'S', [1 0; 0 1; 0 0], ...
                                    'k', k, ...
                                    'nExt', nExt), ...
                'spinw:genmagstr:UCExtNonSuff')
            expected_mag_str = struct('nExt', int32(nExt), ...
                                      'k', k', ...
                                      'F', [1 -1i; 1i 1; 0 0]);
            testCase.verify_obj(swobj.mag_str, expected_mag_str);
        end
        function test_helical_spin_size_incomm_with_epsilon_warns(testCase)
            swobj = copy(testCase.swobj);
            delta = 1e-6;
            k = [(delta+1)/3 0 0];
            nExt = [3 1 1];
            testCase.verifyWarning(...
                @() swobj.genmagstr('mode', 'helical', ...
                                    'S', [1; 0; 0], ...
                                    'k', k, ...
                                    'nExt', nExt, ...
                                    'epsilon', 0.99*delta), ...
                'spinw:genmagstr:UCExtNonSuff')
            expected_mag_str = struct('nExt', int32(nExt), ...
                                      'k', k', ...
                                      'F', [1 -0.5-0.866i -0.5+0.866i; ...
                                            1i 0.866-0.5i -0.866-0.5i; ...
                                             0          0          0]);
            testCase.verify_obj(swobj.mag_str, expected_mag_str, 'rel_tol', 1e-4);
        end
        function test_fourier_too_large_nExt_warns(testCase)
            swobj = copy(testCase.swobj);
            nExt = [6 1 1];
            k = [1/3 0 0];
            testCase.verifyWarning(...
                @() swobj.genmagstr('mode', 'fourier', ...
                                    'S', [1; 1i; 0], ...
                                    'k', k, ...
                                    'nExt', nExt), ...
                'spinw:genmagstr:UCExtOver')
            F_rep = [ 1 -0.5-1i*sqrt(3)/2 -0.5+1i*sqrt(3)/2; ...
                     1i    sqrt(3)/2-0.5i   -sqrt(3)/2-0.5i; ...
                      0                 0                 0];
            expected_mag_str = struct('nExt', int32(nExt), ...
                                      'k', k', ...
                                      'F', cat(2, F_rep, F_rep));
            testCase.verify_obj(swobj.mag_str, expected_mag_str);
        end
        function test_fourier_nExt_wrong_direction_warns(testCase)
            swobj = copy(testCase.swobj);
            nExt = [2 1 2];
            k = [1/2 0 0];
            testCase.verifyWarning(...
                @() swobj.genmagstr('mode', 'fourier', ...
                                    'S', [1; 1i; 0], ...
                                    'k', k, ...
                                    'nExt', nExt), ...
                'spinw:genmagstr:UCExtOver')
            F_rep = [1 -1; 1i -1i; 0 0];
            expected_mag_str = struct('nExt', int32(nExt), ...
                                      'k', k', ...
                                      'F', cat(2, F_rep, F_rep));
            testCase.verify_obj(swobj.mag_str, expected_mag_str);
        end
        function test_helical_any_S_parallel_to_n_warns(testCase)
            swobj_tri = copy(testCase.swobj_tri);
            k = [1/3 0 0];
            testCase.verifyWarning(...
                @() swobj_tri.genmagstr('mode', 'helical', ...
                                        'S', [0; 1; 1], ...
                                        'k', k), ...
                'spinw:genmagstr:SnParallel')
            expected_mag_str = struct( ...
                'nExt', int32([1 1 1]), ...
                'k', k', ...
                'F', [-sqrt(9/8)*1i; sqrt(9/8); sqrt(9/8)]);
            testCase.verify_obj(swobj_tri.mag_str, expected_mag_str);
        end
        function test_helical_S2_norm(testCase, input_norm_output_F)
            swobj = spinw();
            swobj.genlattice('lat_const', [3 8 8], 'angled', [90 90 90]);
            swobj.addatom('r', [0 0 0], 'S', 2);
            k = [1/3 0 0];
            swobj.genmagstr('mode', 'helical', 'S', [1; 0; 0], 'k', k, ...
                            'n', [0 1 0], 'norm', input_norm_output_F{1});
            expected_mag_str = testCase.default_mag_str;
            expected_mag_str.k = k';
            expected_mag_str.F = input_norm_output_F{2};
            testCase.verify_obj(swobj.mag_str, expected_mag_str);
        end
        function test_direct_fm_chain(testCase)
            swobj = copy(testCase.swobj);
            swobj.genmagstr('mode', 'direct', 'k', [0 0 0], ...
                            'S', [0; 1; 0]);
            testCase.verify_obj(swobj.mag_str, testCase.default_mag_str);

        end
        function test_direct_fm_chain_nok(testCase)
            swobj = copy(testCase.swobj);
            swobj.genmagstr('mode', 'direct', 'S', [0; 1; 0]);
            testCase.verify_obj(swobj.mag_str, testCase.default_mag_str);
        end
        function test_direct_multiatom_nExt(testCase)
            swobj = copy(testCase.swobj);
            swobj.addatom('r', [0.5 0.5 0.5], 'S', 2);
            S = [1  0  1 -1; ...
                 1  1  0  0; ...
                 0 -1  0  0];
            nExt = [2 1 1];
            k = [0 1/3 0];
            swobj.genmagstr('mode', 'direct', 'S', S, 'nExt', nExt, 'k', k);
            expected_mag_str = struct('nExt', int32(nExt), ...
                                      'k', k', ...
                                      'F', [sqrt(2)/2        0 1 -2; ...
                                            sqrt(2)/2  sqrt(2) 0  0; ...
                                                    0 -sqrt(2) 0  0]);
            testCase.verify_obj(swobj.mag_str, expected_mag_str);
        end
        function test_direct_multiatom_multik(testCase)
            swobj = copy(testCase.swobj);
            swobj.addatom('r', [0.5 0.5 0.5], 'S', 2);
            S_k = [1 0; 1 1; 0 -1];
            S = cat(3, S_k, S_k);
            k = [0 1/3 0; 1/2 0 0];
            swobj.genmagstr('mode', 'direct', 'S', S, 'k', k);
            F_k = [sqrt(2)/2        0; ...
                   sqrt(2)/2  sqrt(2); ...
                           0 -sqrt(2)];
            expected_mag_str = testCase.default_mag_str;
            expected_mag_str.k = k';
            expected_mag_str.F = cat(3, F_k, F_k);
            testCase.verify_obj(swobj.mag_str, expected_mag_str);
        end
        function test_direct_multik_scalar_nExt(testCase)
            % Test if a scalar is used for nExt it is treated as a
            % tolerance to automatically determine nExt
            swobj = copy(testCase.swobj);
            S_k = [1 0 1 0 1 0; 0 0 1 1 0 0; 0 0 0 1 1 1];
            S = cat(3, S_k, S_k);
            nExt = 0.01;
            k = [0 1/3+0.5*nExt 0; 1/2+0.5*nExt 0 0];
            swobj.genmagstr('mode', 'direct', 'S', S, 'k', k, 'nExt', nExt);
            F_k = [1 0 sqrt(2)/2         0 sqrt(2)/2 0; ...
                   0 0 sqrt(2)/2 sqrt(2)/2         0 0; ...
                   0 0         0 sqrt(2)/2 sqrt(2)/2 1];
            expected_mag_str = struct('nExt', int32([2 3 1]), ...
                                      'k', k', ...
                                      'F', cat(3, F_k, F_k));
            testCase.verify_obj(swobj.mag_str, expected_mag_str);
        end
        function test_helical_tri(testCase)
            swobj_tri = copy(testCase.swobj_tri);
            k = [1/3 1/3 0];
            swobj_tri.genmagstr('mode', 'helical', 'S', [1; 0; 0], ...
                                'k', k);
            expected_mag_str = testCase.default_mag_str;
            expected_mag_str.k = k';
            expected_mag_str.F = [1.5; 1.5i; 0];
            testCase.verify_obj(swobj_tri.mag_str, expected_mag_str);
        end
        function test_helical_tri_n(testCase)
            swobj_tri = copy(testCase.swobj_tri);
            k = [1/3 1/3 0];
            swobj_tri.genmagstr('mode', 'helical', 'S', [1; 0; 0], ...
                                'k', k, 'n', [0 1 0]);
            expected_mag_str = testCase.default_mag_str;
            expected_mag_str.k = k';
            expected_mag_str.F = [1.5; 0; -1.5i];
            testCase.verify_obj(swobj_tri.mag_str, expected_mag_str);
        end
        function test_helical_tri_lu_unit(testCase)
            swobj_tri = copy(testCase.swobj_tri);
            swobj_tri.genmagstr('mode', 'helical', 'S', [0; 1; 0], ...
                                'k', [1/3 1/3 0], 'unit', 'lu');
            expected_mag_str = testCase.default_mag_str;
            expected_mag_str.k = [1/3; 1/3; 0];
            expected_mag_str.F = [-0.75-1.299038105676658i; ...
                                  1.299038105676658-0.75i; ...
                                  0];
            testCase.verify_obj(swobj_tri.mag_str, expected_mag_str);
        end
        function test_fourier_tri(testCase)
            swobj_tri = copy(testCase.swobj_tri);
            swobj_tri.genmagstr('mode', 'fourier', 'S', [1; 0; 0], ...
                                'k', [1/3 1/3 0]);
            expected_mag_str = testCase.default_mag_str;
            expected_mag_str.k = [1/3; 1/3; 0];
            expected_mag_str.F = [1.5; 0; 0];
            testCase.verify_obj(swobj_tri.mag_str, expected_mag_str);
        end
        function test_helical_multiatom_nExt_1spin(testCase)
            % Test where only 1 spin is provided
            swobj = copy(testCase.swobj);
            swobj.addatom('r', [0.5 0.5 0.0], 'S', 2);
            k = [0 0 1/2];
            nExt = int32([1 1 2]);
            swobj.genmagstr('mode', 'helical', 'S', [1; 0; 0], ...
                            'nExt', nExt, 'k', k);
            expected_mag_str = struct(...
                'k', k', ...
                'nExt', nExt, ...
                'F', [1 2 -1 -2; 1i 2i -1i -2i; 0 0 0 0]);
             testCase.verify_obj(swobj.mag_str, expected_mag_str);
        end
        function test_helical_multiatom_nExt_1spin_r0(testCase)
            swobj = spinw();
            swobj.genlattice('lat_const', [3 8 8], 'angled', [90 90 90]);
            % Need to have nonzero r for first atom for r0 to have an effect
            swobj.addatom('r', [0.5 0.5 0.5],'S', 1);
            swobj.addatom('r', [0 0 0], 'S', 2);
            k = [0 0 1/2];
            nExt = int32([1 1 2]);
            swobj.genmagstr('mode', 'helical', 'S', [1; 0; 0], ...
                            'nExt', nExt, 'k', k, 'r0', false);
            expected_mag_str = struct(...
                'k', k', ...
                'nExt', nExt, ...
                'F', [1 2i -1 -2i; 1i -2 -1i 2; 0 0 0 0]);
             testCase.verify_obj(swobj.mag_str, expected_mag_str);
        end
        function test_fourier_multiatom_nExt_nMagAtom_spins(testCase)
            % Test where there are the same number of spins provided as in
            % the unit cell. Note result is the same as helical with
            % complex spins or S=[0 1; 1 0; 0 0]
            swobj = copy(testCase.swobj);
            swobj.addatom('r', [0.5 0.5 0.0], 'S', 2);
            k = [0 1/2 0];
            nExt = int32([1 2 1]);
            swobj.genmagstr('mode', 'fourier', 'S', [-1i 1; 1 1i; 0 0], ...
                            'nExt', nExt, 'k', k);
            expected_mag_str = struct( ...
                'k', k', ...
                'nExt', nExt, ...
                'F', [-1i 2 1i -2; 1 2i -1 -2i; 0 0 0 0]);
             testCase.verify_obj(swobj.mag_str, expected_mag_str);
        end
        function test_helical_multiatom_nExt_nMagAtom_spins(testCase)
            % Test where there are the same number of spins provided as in
            % the unit cell
            swobj = copy(testCase.swobj);
            swobj.addatom('r', [0.5 0.5 0.0], 'S', 2);
            k = [0 1/2 0];
            nExt = int32([1 2 1]);
            swobj.genmagstr('mode', 'helical', 'S', [0 1; 1 0; 0 0], ...
                            'nExt', nExt, 'k', k);
            expected_mag_str = struct( ...
                'k', k', ...
                'nExt', nExt, ...
                'F', [-1i 2 1i -2; 1 2i -1 -2i; 0 0 0 0]);
             testCase.verify_obj(swobj.mag_str, expected_mag_str);
        end
        function test_helical_multiatom_nExt_nMagExt_spins(testCase)
            % Test where there are the same number of spins provided as in
            % the supercell
            swobj = copy(testCase.swobj);
            swobj.addatom('r', [0.5 0.5 0.0], 'S', 2);
            k = [0 1/2 0];
            nExt = int32([1 2 1]);
            testCase.verifyWarning(...
                @() swobj.genmagstr( ...
                    'mode', 'helical', ...
                    'S', [0 1 0 -1; 1 0 0 0; 0 0 1 0], ...
                    'nExt', nExt, 'k', k), ...
                'spinw:genmagstr:SnParallel');
            expected_mag_str = struct(...
                'k', k', ...
                'nExt', nExt, ...
                'F', [-1i 2 0 -2; 1 2i 0 -2i; 0 0 1 0]);
             testCase.verify_obj(swobj.mag_str, expected_mag_str);
        end
        function test_helical_multiatom_multik_multin(testCase)
            swobj = copy(testCase.swobj);
            swobj.addatom('r', [0.5 0.5 0.5], 'S', 2);
            S = cat(3, [1; 0; 0], [0; 0; 1]);
            k = [0 1/3 0; 1/2 0 0];
            n = [0 0 1; 0 1 0];
            % Ensure warning is not emitted as there are no S parallel to
            % n within a single k
            testCase.verifyWarningFree(...
                @() swobj.genmagstr('mode', 'helical', ...
                                    'S', S, ...
                                    'k', k, ...
                                    'n', n), ...
                'spinw:genmagstr:SnParallel')
            expected_mag_str = testCase.default_mag_str;
            expected_mag_str.k = k';
            expected_mag_str.F = cat(3, ...
                                     [1 1-sqrt(3)*1i; 1i sqrt(3)+1i; 0   0], ...
                                     [1i           2;  0          0; 1 -2i]);
            testCase.verify_obj(swobj.mag_str, expected_mag_str);
        end
        function test_random_structure(testCase)
            swobj = copy(testCase.swobj);
            swobj.genmagstr('mode','random');
            mag_str1 = swobj.mag_str;
            swobj.genmagstr('mode','random');
            mag_str2 = swobj.mag_str;
            % Check structure is random each time - F is different
            testCase.verifyNotEqual(mag_str1.F, mag_str2.F);
            testCase.verifyEqual(size(mag_str1.F), size(mag_str2.F));
            % Check F size and magnitude
            testCase.verifySize(mag_str1.F, [3 1]);
            testCase.verify_val(vecnorm(real(swobj.mag_str.F), 2), 1);
            % Check imaginary component of F is perpendicular to default n
            testCase.verify_val(dot(imag(swobj.mag_str.F), [0 0 1]), 0);
            % Check other fields
            expected_mag_str = testCase.default_mag_str;
            testCase.verify_obj(rmfield(mag_str1, 'F'), ...
                                rmfield(expected_mag_str, 'F'));
        end
        function test_random_structure_k_and_n(testCase)
            swobj = copy(testCase.swobj);
            k = [0; 0; 1/4];
            n = [1 1 0];
            swobj.genmagstr('mode','random', 'k', k', 'n', n);
            mag_str1 = swobj.mag_str;
            % Check F size and magnitude
            testCase.verifySize(mag_str1.F, [3 1]);
            testCase.verify_val(vecnorm(real(swobj.mag_str.F), 2), 1);
            % Check imaginary component of F is perpendicular to n
            testCase.verify_val(dot(imag(swobj.mag_str.F), n), 0);
            % Check other fields
            expected_mag_str = testCase.default_mag_str;
            expected_mag_str.k = k;
            testCase.verify_obj(rmfield(mag_str1, 'F'),...
                                rmfield(expected_mag_str, 'F'));
        end
        function test_random_structure_multiatom_and_nExt(testCase)
            swobj = copy(testCase.swobj);
            swobj.addatom('r', [0.5 0.5 0], 'S', 2);
            nExt = int32([2 2 2]);
            swobj.genmagstr('mode', 'random', 'nExt', nExt);
            mag_str1 = swobj.mag_str;
            swobj.genmagstr('mode', 'random', 'nExt', nExt);
            mag_str2 = swobj.mag_str;
            % Check structure is random each time - F is different
            testCase.verifyNotEqual(mag_str1.F, mag_str2.F);
            % Check F size and magnitude
            testCase.verifySize(mag_str1.F, [3 16]);
            testCase.verify_val( ...
                vecnorm(real(swobj.mag_str.F(:, 1:2:end)), 2), ones(1, 8));
            testCase.verify_val( ...
                vecnorm(real(swobj.mag_str.F(:, 2:2:end)), 2), 2*ones(1, 8));
            % Check imaginary component of F is perpendicular to default n
            testCase.verifyEqual( ...
                dot(imag(swobj.mag_str.F), repmat([0; 0; 1], 1, 16)), zeros(1, 16));
            % Check other fields
            expected_mag_str = testCase.default_mag_str;
            expected_mag_str.nExt = nExt;
            testCase.verify_obj(rmfield(mag_str1, 'F'), ...
                                rmfield(expected_mag_str, 'F'));
        end
        function test_tile_existing_struct_extend_cell(testCase)
            % Test that tile and increasing nExt will correctly tile
            % initialised structure
            swobj = copy(testCase.swobj);
            swobj.addatom('r', [0.5 0.5 0], 'S', 1);
            nExt = int32([1 2 1]);
            % Also test if we input 'k' it is set to 0 in final struct
            swobj.genmagstr('mode', 'direct', 'S', [1 0; 0 1; 0 0], 'k', [1/2 0 0]);
            swobj.genmagstr('mode', 'tile', 'nExt', nExt);
            expected_mag_str = testCase.default_mag_str;
            expected_mag_str.nExt = nExt;
            expected_mag_str.F = [1 0 1 0; 0 1 0 1; 0 0 0 0];
            testCase.verify_obj(swobj.mag_str, expected_mag_str);
        end
        function test_tile_existing_struct_same_size(testCase)
            % Test that tile with nExt same as initialised structure
            % does nothing
            swobj = copy(testCase.swobj);
            swobj.addatom('r', [0.5 0.5 0], 'S', 1);
            nExt = int32([1 2 1]);
            S = [1 0 0 -1; 0 1 0 0; 0 0 1 0];
            swobj.genmagstr('mode', 'direct', 'S', S, 'nExt', nExt);
            swobj.genmagstr('mode', 'tile', 'nExt', nExt);
            expected_mag_str = testCase.default_mag_str;
            expected_mag_str.nExt = nExt;
            expected_mag_str.F = S;
            testCase.verify_obj(swobj.mag_str, expected_mag_str);
        end
        function test_tile_input_S_extend_cell(testCase)
            % Test that tile and input S less than nExt will correctly tile
            % input S
            swobj = copy(testCase.swobj);
            swobj.addatom('r', [0.5 0.5 0], 'S', 1);
            nExt = int32([3 1 1]);
            swobj.genmagstr('mode', 'tile', 'nExt', nExt, ...
                            'S', [1 0; 0 1; 0 0]);
            expected_mag_str = testCase.default_mag_str;
            expected_mag_str.nExt = nExt;
            expected_mag_str.F = [1 0 1 0 1 0; 0 1 0 1 0 1; 0 0 0 0 0 0];
            testCase.verify_obj(swobj.mag_str, expected_mag_str);
        end
        function test_tile_multik(testCase)
            % Test that S is summed over third dimension with tile, and k
            % is not needed (is this the behaviour we want?)
            swobj = copy(testCase.swobj);
            swobj.addatom('r', [0.5 0.5 0], 'S', 1);
            S = cat(3, [1 0; 0 1; 0 0], [0 1; 0 0; 1 0]);
            swobj.genmagstr('mode', 'tile', ...
                            'S', S);
            expected_mag_str = testCase.default_mag_str;
            expected_mag_str.F = sqrt(2)/2*[1 1; 0 1; 1 0];
            testCase.verify_obj(swobj.mag_str, expected_mag_str);
        end
        function test_tile_multik_provided_k_set_to_zero(testCase)
            % Test that S is summed over third dimension with tile, and if
            % k is provided, it is set to zero anyway
            swobj = copy(testCase.swobj);
            swobj.addatom('r', [0.5 0.5 0], 'S', 1);
            S = cat(3, [1 0; 0 1; 0 0], [0 1; 0 0; 1 0]);
            k = 0.5*ones(size(S, 3), 3);
            swobj.genmagstr('mode', 'tile', ...
                              'S', S, 'k', k);
            expected_mag_str = testCase.default_mag_str;
            expected_mag_str.F = sqrt(2)/2*[1 1; 0 1; 1 0];
            testCase.verify_obj(swobj.mag_str, expected_mag_str);
        end
        function test_extend_mode_input_S_extend_cell_and_warns(testCase)
            % Test undocumented 'extend' mode does same as tile
            % Test that tile and input S less than nExt will correctly tile
            % input S
            swobj = copy(testCase.swobj);
            swobj.addatom('r', [0.5 0.5 0], 'S', 1);
            nExt = int32([3 1 1]);
            testCase.verifyWarning(...
                @() swobj.genmagstr('mode', 'extend', 'nExt', nExt, ...
                                    'S', [1 0; 0 1; 0 0]), ...
                'spinw:genmagstr:DeprecationWarning');
            expected_mag_str = testCase.default_mag_str;
            expected_mag_str.nExt = nExt;
            expected_mag_str.F = [1 0 1 0 1 0; 0 1 0 1 0 1; 0 0 0 0 0 0];
            testCase.verify_obj(swobj.mag_str, expected_mag_str);
        end
        function test_rotate_phi(testCase, rotate_phi_inputs)
            swobj = copy(testCase.swobj);
            k = [1/2 0 0];
            % Need to initialise structure before rotating it
            swobj.genmagstr('mode', 'direct', 'S', [1; 0; 0], 'k', k);
            swobj.genmagstr('mode', 'rotate', rotate_phi_inputs{:});
            expected_mag_str = testCase.default_mag_str;
            expected_mag_str.F = sqrt(2)/2*[1; 1; 0];
            expected_mag_str.k = k';
            testCase.verify_obj(swobj.mag_str, expected_mag_str);
        end
        function test_rotate_multiatom_n(testCase)
            swobj = copy(testCase.swobj);
            swobj.addatom('r', [0.5 0.5 0], 'S', 1);
            k = [1/2 0 0];
            swobj.genmagstr('mode', 'direct', 'S', [1 0; 0 1; 0 0], 'k', k);
            swobj.genmagstr('mode', 'rotate', 'phi', pi/2, 'n', [1 1 0]);
            expected_mag_str = testCase.default_mag_str;
            expected_mag_str.F = [0.5 0.5; 0.5 0.5; -sqrt(2)/2 sqrt(2)/2];
            expected_mag_str.k = k';
            testCase.verify_obj(swobj.mag_str, expected_mag_str);
        end
        function test_rotate_no_phi_collinear(testCase)
            swobj = copy(testCase.swobj);
            swobj.addatom('r', [0.5 0.5 0], 'S', 1);
            swobj.genmagstr('mode', 'direct', 'S', [1 -1; 0 0; 0 0]);
            swobj.genmagstr('mode', 'rotate', 'n', [0 1 0]);
            expected_mag_str = testCase.default_mag_str;
            expected_mag_str.F = [0 0; 1 -1; 0 0];
            testCase.verify_obj(swobj.mag_str, expected_mag_str);
        end
        function test_rotate_no_phi_coplanar(testCase)
            swobj = copy(testCase.swobj);
            swobj.addatom('r', [0.5 0.5 0], 'S', 1);
            swobj.genmagstr('mode', 'direct', 'S', [1 0; 0 1; 0 0]);
            swobj.genmagstr('mode', 'rotate', 'n', [0 1 0]);
            expected_mag_str = testCase.default_mag_str;
            expected_mag_str.F = [1 0; 0 0; 0 -1];
            testCase.verify_obj(swobj.mag_str, expected_mag_str);
        end
        function test_rotate_no_phi_incomm(testCase)
            swobj_tri = copy(testCase.swobj_tri);
            k = [1/3 1/3 0];
            swobj_tri.genmagstr('mode', 'helical', 'S', [1; 0; 0], ...
                                'k', k);
            swobj_tri.genmagstr('mode', 'rotate', 'n', [1 0 0]);
            expected_mag_str = testCase.default_mag_str;
            expected_mag_str.F = [0; 1.5i; -1.5];
            expected_mag_str.k = k';
            testCase.verify_obj(swobj_tri.mag_str, expected_mag_str);
        end
        function test_func_multiatom_default(testCase)
            swobj = copy(testCase.swobj);
            swobj.addatom('r', [0.5 0.5 0], 'S', 1);
            k = [1/3 0 0];
            x0 = [pi/2 -pi/4 0 0 k pi pi/2];
            swobj.genmagstr('mode', 'func', 'x0', x0);
            expected_mag_str = testCase.default_mag_str;
            expected_mag_str.F = [sqrt(2)/2*(1-i) 0; -sqrt(2)/2*(1+i) 0; 0 1];
            expected_mag_str.k = k';
            testCase.verify_obj(swobj.mag_str, expected_mag_str);
        end
        function test_func_custom(testCase)
             function [S, k, n] = func(S0, x0)
                 S = [-S0; 0; 0];
                 k = x0;
                 n = x0;
             end
            swobj = copy(testCase.swobj);
            x0 = [1/3 0 0];
            swobj.genmagstr('mode', 'func', 'func', @func, 'x0', x0);
            expected_mag_str = testCase.default_mag_str;
            expected_mag_str.F = [-1; 0; 0];
            expected_mag_str.k = x0';
            testCase.verify_obj(swobj.mag_str, expected_mag_str);
        end
    end
    methods (Test, TestTags = {'Symbolic'})
        function test_func_custom_symbolic(testCase)
            function [S, k, n] = func(S0, x0)
                S = [-S0; 0; 0];
                k = x0;
                n = x0;
            end
            swobj = copy(testCase.swobj);
            swobj.symbolic(true);
            x0 = [1/3 0 0];
            swobj.genmagstr('mode', 'func', 'func', @func, 'x0', x0);
            expected_mag_str = testCase.default_mag_str;
            expected_mag_str.F = sym([-1; 0; 0]);
            expected_mag_str.k = sym(x0');
            testCase.verify_obj(swobj.mag_str, expected_mag_str);
        end
    end
end