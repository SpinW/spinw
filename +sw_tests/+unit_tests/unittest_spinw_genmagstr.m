classdef unittest_spinw_genmagstr < sw_tests.unit_tests.unittest_super

    properties
        swobj = [];
        swobj_tri = [];
        default_mag_str = struct('nExt', ones(1, 3, 'int32'), ...
                                 'k', [0; 0; 0], ...
                                 'F', [0; 1; 0]);
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
        function test_no_magnetic_atoms_raises_error(testCase)
            swobj = spinw;
            testCase.verifyError(...
                @() swobj.genmagstr(), ...
                'spinw:genmagstr:NoMagAtom')
        end
        function test_invalid_mode_raises_error(testCase)
            testCase.verifyError(...
                @() testCase.swobj.genmagstr('mode', 'something'), ...
                'spinw:genmagstr:WrongMode')
        end
        function test_zero_nExt_raises_error(testCase)
            testCase.verifyError(...
                @() testCase.swobj.genmagstr('nExt', [0 1 1]), ...
                'spinw:genmagstr:WrongInput')
        end
        function test_direct_wrong_spin_size_raises_error(testCase)
            testCase.verifyError(...
                @() testCase.swobj.genmagstr( ...
                    'mode', 'direct', 'S', [0; 1; 0], 'nExt', [2 1 1]), ...
                'spinw:genmagstr:WrongSpinSize')
        end
        function test_helical_wrong_spin_size_warns(testCase)
            testCase.verifyWarning(...
                @() testCase.swobj.genmagstr('mode', 'helical', ...
                                             'S', [1 0; 0 1; 0 0], ...
                                             'k', [1/3 0 0], ...
                                             'nExt', [2 1 1]), ...
                'spinw:genmagstr:UCExtNonSuff')
        end
        function test_direct_fm_chain(testCase)
            swobj = copy(testCase.swobj);
            swobj.genmagstr('mode', 'direct', 'k', [0 0 0], ...
                            'S', [0; 1; 0]);
            testCase.verify_obj(testCase.default_mag_str, swobj.mag_str);

        end
        function test_direct_fm_chain_nok(testCase)
            swobj = copy(testCase.swobj);
            swobj.genmagstr('mode', 'direct', 'S', [0; 1; 0]);
            testCase.verify_obj(testCase.default_mag_str, swobj.mag_str);
        end
        function test_direct_multiatom_nExt(testCase)
            swobj = copy(testCase.swobj);
            swobj.addatom('r', [0.5 0.5 0.5], 'S', 2);
            S = [1  0  1 -1; ...
                 1  1  0  0; ...
                 0 -1  0  0];
            nExt = [int32(2) int32(1) int32(1)];
            k = [0 1/3 0];
            swobj.genmagstr('mode', 'direct', 'S', S, 'nExt', nExt, 'k', k);
            expected_mag_str = struct('nExt', nExt, ...
                                      'k', k', ...
                                      'F', [sqrt(2)/2        0 1 -2; ...
                                            sqrt(2)/2  sqrt(2) 0  0; ...
                                                    0 -sqrt(2) 0  0]);
            testCase.verify_obj(expected_mag_str, swobj.mag_str);
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
            testCase.verify_obj(expected_mag_str, swobj.mag_str);
        end
        function test_helical_tri(testCase)
            swobj_tri = copy(testCase.swobj_tri);
            swobj_tri.genmagstr('mode', 'helical', 'S', [1; 0; 0], ...
                                'k', [1/3 1/3 0]);
            expected_mag_str = testCase.default_mag_str;
            expected_mag_str.k = [1/3; 1/3; 0];
            expected_mag_str.F = [1.5; 1.5i; 0];
            testCase.verify_obj(expected_mag_str, swobj_tri.mag_str);
        end
        function test_fourier_tri(testCase)
            swobj_tri = copy(testCase.swobj_tri);
            swobj_tri.genmagstr('mode', 'fourier', 'S', [1; 0; 0], ...
                                'k', [1/3 1/3 0]);
            expected_mag_str = testCase.default_mag_str;
            expected_mag_str.k = [1/3; 1/3; 0];
            expected_mag_str.F = [1.5; 0; 0];
            testCase.verify_obj(expected_mag_str, swobj_tri.mag_str);
        end
        function test_helical_multiatom_nExt_1spin(testCase)
            % Test where only 1 spin is provided
            swobj = copy(testCase.swobj);
            swobj.addatom('r', [0.5 0.5 0.0], 'S', 2);
            k = [0 0 1/2];
            nExt = [int32(1) int32(1) int32(2)];
            swobj.genmagstr('mode', 'helical', 'S', [1; 0; 0], ...
                            'nExt', nExt, 'k', k);
            expected_mag_str = struct(...
                'k', k', ...
                'nExt', nExt, ...
                'F', [1 2 -1 -2; 1i 2i -1i -2i; 0 0 0 0]);
             testCase.verify_obj(expected_mag_str, swobj.mag_str);
        end
        function test_helical_multiatom_nExt_1spin_r0(testCase)
            swobj = spinw();
            swobj.genlattice('lat_const', [3 8 8], 'angled', [90 90 90]);
            % Need to have nonzero r for first atom for r0 to have an effect
            swobj.addatom('r', [0.5 0.5 0.5],'S', 1);
            swobj.addatom('r', [0 0 0], 'S', 2);
            k = [0 0 1/2];
            nExt = [int32(1) int32(1) int32(2)];
            swobj.genmagstr('mode', 'helical', 'S', [1; 0; 0], ...
                            'nExt', nExt, 'k', k, 'r0', false);
            expected_mag_str = struct(...
                'k', k', ...
                'nExt', nExt, ...
                'F', [1 2i -1 -2i; 1i -2 -1i 2; 0 0 0 0]);
             testCase.verify_obj(expected_mag_str, swobj.mag_str);
        end
        function test_fourier_multiatom_nExt_nMagAtom_spins(testCase)
            % Test where there are the same number of spins provided as in
            % the unit cell. Note result is the same as helical with
            % complex spins or S=[0 1; 1 0; 0 0]
            swobj = copy(testCase.swobj);
            swobj.addatom('r', [0.5 0.5 0.0], 'S', 2);
            k = [0 1/2 0];
            nExt = [int32(1) int32(2) int32(1)];
            swobj.genmagstr('mode', 'fourier', 'S', [-1i 1; 1 1i; 0 0], ...
                            'nExt', nExt, 'k', k);
            expected_mag_str = struct( ...
                'k', k', ...
                'nExt', nExt, ...
                'F', [-1i 2 1i -2; 1 2i -1 -2i; 0 0 0 0]);
             testCase.verify_obj(expected_mag_str, swobj.mag_str);
        end
        function test_helical_multiatom_nExt_nMagAtom_spins(testCase)
            % Test where there are the same number of spins provided as in
            % the unit cell
            swobj = copy(testCase.swobj);
            swobj.addatom('r', [0.5 0.5 0.0], 'S', 2);
            k = [0 1/2 0];
            nExt = [int32(1) int32(2) int32(1)];
            swobj.genmagstr('mode', 'helical', 'S', [0 1; 1 0; 0 0], ...
                            'nExt', nExt, 'k', k);
            expected_mag_str = struct( ...
                'k', k', ...
                'nExt', nExt, ...
                'F', [-1i 2 1i -2; 1 2i -1 -2i; 0 0 0 0]);
             testCase.verify_obj(expected_mag_str, swobj.mag_str);
        end
        function test_helical_multiatom_nExt_nMagExt_spins(testCase)
            % Test where there are the same number of spins provided as in
            % the supercell
            swobj = copy(testCase.swobj);
            swobj.addatom('r', [0.5 0.5 0.0], 'S', 2);
            k = [0 1/2 0];
            nExt = [int32(1) int32(2) int32(1)];
            swobj.genmagstr('mode', 'helical', ...
                            'S', [0 1 0 -1; 1 0 0 0; 0 0 1 0], ...
                            'nExt', nExt, 'k', k);
            expected_mag_str = struct(...
                'k', k', ...
                'nExt', nExt, ...
                'F', [-1i 2 0 -2; 1 2i 0 -2i; 0 0 1 0]);
             testCase.verify_obj(expected_mag_str, swobj.mag_str);
        end
        function test_helical_multiatom_multik_multin(testCase)
            swobj = copy(testCase.swobj);
            swobj.addatom('r', [0.5 0.5 0.5], 'S', 2);
            S_k = [1; 0; 0];
            S = cat(3, S_k, S_k);
            k = [0 1/3 0; 1/2 0 0];
            n = [0 0 1; 0 1 0];
            swobj.genmagstr('mode', 'helical', 'S', S, 'k', k, 'n', n);
            F_k = [sqrt(2)/2        0; ...
                   sqrt(2)/2  sqrt(2); ...
                           0 -sqrt(2)];
            expected_mag_str = testCase.default_mag_str;
            expected_mag_str.k = k';
            expected_mag_str.F = cat(3, ...
                                     [1 1-sqrt(3)*1i; 1i sqrt(3)+1i;   0  0], ...
                                     [1         -2i;  0           0; -1i -2]);

            testCase.verify_obj(expected_mag_str, swobj.mag_str);
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
            % Check S has correct size and magnitude
            testCase.verify_val(norm(swobj.magstr.S, 2), 1);
            % Check other fields
            expected_mag_str = testCase.default_mag_str;
            testCase.verify_obj(rmfield(expected_mag_str, 'F'), ...
                                rmfield(mag_str1, 'F'));
        end
        function test_random_structure_k(testCase)
            swobj = copy(testCase.swobj);
            k = [0; 0; 1/4];
            swobj.genmagstr('mode','random', 'k', k');
            mag_str1 = swobj.mag_str;
            % Check F size
            testCase.verifySize(mag_str1.F, [3 1]);
            % Check S has correct size and magnitude
            testCase.verify_val(vecnorm(swobj.magstr.S, 2), 1);
            % Check other fields
            expected_mag_str = testCase.default_mag_str;
            expected_mag_str.k = k;
            testCase.verify_obj(rmfield(expected_mag_str, 'F'), ...
                                rmfield(mag_str1, 'F'));
        end
        function test_random_structure_multiatom_and_nExt(testCase)
            swobj = copy(testCase.swobj);
            swobj.addatom('r', [0.5 0.5 0], 'S', 2);
            nExt = 2*ones(1, 3, 'int32');
            swobj.genmagstr('mode', 'random', 'nExt', nExt);
            mag_str1 = swobj.mag_str;
            swobj.genmagstr('mode', 'random', 'nExt', nExt);
            mag_str2 = swobj.mag_str;
            % Check structure is random each time - F is different
            testCase.verifyNotEqual(mag_str1.F, mag_str2.F);
            testCase.verifySize(mag_str1.F, [3 16]);
            % Check S has correct magnitudes
            testCase.verify_val(vecnorm(swobj.magstr.S(:, 1:2:end), 2), ...
                                ones(1, 8));
            testCase.verify_val(vecnorm(swobj.magstr.S(:, 2:2:end), 2), ...
                                2*ones(1, 8));
            % Check other fields
            expected_mag_str = testCase.default_mag_str;
            expected_mag_str.nExt = nExt;
            testCase.verify_obj(rmfield(expected_mag_str, 'F'), ...
                                rmfield(mag_str1, 'F'));
        end
        function test_tile_no_init_creates_random(testCase)
            % Test tile without initialised structure
            % (i.e. size(S,2) < nMagAtom) creates a random structure
            swobj1 = copy(testCase.swobj);
            swobj1.addatom('r', [0.5 0.5 0], 'S', 1);
            % Need to use clean obj or the stored structure is used
            swobj2 = copy(swobj1);
            nExt = [int32(2) int32(1) int32(1)];
            swobj1.genmagstr('mode', 'tile', 'nExt', nExt, 'S', [1; 0; 0]);
            mag_str1 = swobj1.mag_str;
            swobj2.genmagstr('mode', 'tile', 'nExt', nExt, 'S', [1; 0; 0]);
            mag_str2 = swobj2.mag_str;
            % Check structure is random each time - F is different
            testCase.verifyNotEqual(mag_str1.F, mag_str2.F);
            testCase.verifySize(mag_str1.F, [3 4]);
            % Check S has correct magnitudes
            testCase.verify_val(vecnorm(swobj1.magstr.S, 2), ones(1, 4));
            % Check other fields
            expected_mag_str = testCase.default_mag_str;
            expected_mag_str.nExt = nExt;
            testCase.verify_obj(rmfield(expected_mag_str, 'F'), ...
                                rmfield(mag_str1, 'F'));
        end
        function test_tile_existing_struct_extend_cell(testCase)
            % Test that tile and increasing nExt will correctly tile
            % initialised structure
            swobj = copy(testCase.swobj);
            swobj.addatom('r', [0.5 0.5 0], 'S', 1);
            nExt = [int32(1) int32(2) int32(1)];
            % Also test if we input 'k' it is set to 0 in final struct
            swobj.genmagstr('mode', 'direct', 'S', [1 0; 0 1; 0 0], 'k', [1/2 0 0]);
            swobj.genmagstr('mode', 'tile', 'nExt', nExt);
            expected_mag_str = testCase.default_mag_str;
            expected_mag_str.nExt = nExt;
            expected_mag_str.F = [1 0 1 0; 0 1 0 1; 0 0 0 0];
            testCase.verify_obj(expected_mag_str, swobj.mag_str);
        end
        function test_tile_existing_struct_same_size(testCase)
            % Test that tile with nExt same as initialised structure
            % does nothing
            swobj = copy(testCase.swobj);
            swobj.addatom('r', [0.5 0.5 0], 'S', 1);
            nExt = [int32(1) int32(2) int32(1)];
            S = [1 0 0 -1; 0 1 0 0; 0 0 1 0];
            swobj.genmagstr('mode', 'direct', 'S', S, 'nExt', nExt);
            swobj.genmagstr('mode', 'tile', 'nExt', nExt);
            expected_mag_str = testCase.default_mag_str;
            expected_mag_str.nExt = nExt;
            expected_mag_str.F = S;
            testCase.verify_obj(expected_mag_str, swobj.mag_str);
        end
        function test_tile_input_S_extend_cell(testCase)
            % Test that tile and input S less than nExt will correctly tile
            % input S
            swobj = copy(testCase.swobj);
            swobj.addatom('r', [0.5 0.5 0], 'S', 1);
            nExt = [int32(3) int32(1) int32(1)];
            swobj.genmagstr('mode', 'tile', 'nExt', nExt, ...
                            'S', [1 0; 0 1; 0 0]);
            expected_mag_str = testCase.default_mag_str;
            expected_mag_str.nExt = nExt;
            expected_mag_str.F = [1 0 1 0 1 0; 0 1 0 1 0 1; 0 0 0 0 0 0];
            testCase.verify_obj(expected_mag_str, swobj.mag_str);
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
            testCase.verify_obj(expected_mag_str, swobj.mag_str);
        end
        function test_rotate_phi(testCase)
            swobj = copy(testCase.swobj);
            k = [1/2 0 0];
            % Need to initialise structure before rotating it
            swobj.genmagstr('mode', 'direct', 'S', [1; 0; 0], 'k', k);
            swobj.genmagstr('mode', 'rotate', 'phi', pi/4);
            expected_mag_str = testCase.default_mag_str;
            expected_mag_str.F = sqrt(2)/2*[1; 1; 0];
            expected_mag_str.k = k';
            testCase.verify_obj(expected_mag_str, swobj.mag_str);
        end
        function test_rotate_phid(testCase)
            swobj = copy(testCase.swobj);
            k = [1/2 0 0];
            swobj.genmagstr('mode', 'direct', 'S', [1; 0; 0], 'k', k);
            swobj.genmagstr('mode', 'rotate', 'phid', 45);
            expected_mag_str = testCase.default_mag_str;
            expected_mag_str.F = sqrt(2)/2*[1; 1; 0];
            expected_mag_str.k = k';
            testCase.verify_obj(expected_mag_str, swobj.mag_str);
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
            testCase.verify_obj(expected_mag_str, swobj.mag_str);
        end
        function test_rotate_no_phi_collinear(testCase)
            swobj = copy(testCase.swobj);
            swobj.addatom('r', [0.5 0.5 0], 'S', 1);
            swobj.genmagstr('mode', 'direct', 'S', [1 -1; 0 0; 0 0]);
            swobj.genmagstr('mode', 'rotate', 'n', [0 1 0]);
            expected_mag_str = testCase.default_mag_str;
            expected_mag_str.F = [0 0; 1 -1; 0 0];
            testCase.verify_obj(expected_mag_str, swobj.mag_str);
        end
        function test_rotate_no_phi_coplanar(testCase)
            swobj = copy(testCase.swobj);
            swobj.addatom('r', [0.5 0.5 0], 'S', 1);
            swobj.genmagstr('mode', 'direct', 'S', [1 0; 0 1; 0 0]);
            swobj.genmagstr('mode', 'rotate', 'n', [0 1 0]);
            expected_mag_str = testCase.default_mag_str;
            expected_mag_str.F = [1 0; 0 0; 0 -1];
            testCase.verify_obj(expected_mag_str, swobj.mag_str);
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
            testCase.verify_obj(expected_mag_str, swobj_tri.mag_str);
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
            testCase.verify_obj(expected_mag_str, swobj.mag_str);
        end
    end

end