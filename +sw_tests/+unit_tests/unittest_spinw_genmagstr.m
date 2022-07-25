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
            testCase.swobj.addatom('r', [0 0 0],'S', 1, 'label', 'MNi2');
        end
        function setup_tri_model(testCase)
            % Create a simple triangular lattice model
            testCase.swobj_tri = spinw;
            testCase.swobj_tri.genlattice('lat_const', [4 4 6], 'angled', [90 90 120]);
            testCase.swobj_tri.addatom('r', [0 0 0], 'S', 3/2, 'label', 'MCr3');
        end
    end

    methods (Test)
        function test_no_magnetic_atoms_raises_error(testCase)
            swobj = spinw;
            testCase.verifyError(...
                @() swobj.genmagstr(), ...
                'spinw:genmagstr:NoMagAtom')
        end
        function test_zero_nExt_raises_error(testCase)
            testCase.verifyError(...
                @() testCase.swobj.genmagstr('nExt', [0 1 1]), ...
                'spinw:genmagstr:WrongInput')
        end
        function test_direct_fm_chain(testCase)
            swobj = copy(testCase.swobj);
            swobj.genmagstr('mode', 'direct', 'k', [0 0 0], ...
                            'n', [1 0 0], 'S', [0; 1; 0]);
            testCase.verify_obj(testCase.default_mag_str, swobj.mag_str);

        end
        function test_direct_fm_chain_nok(testCase)
            swobj = copy(testCase.swobj);
            swobj.genmagstr('mode', 'direct', 'n', [1 0 0], 'S', [0; 1; 0]);
            testCase.verify_obj(testCase.default_mag_str, swobj.mag_str);
        end
        function test_helical_tri(testCase)
            swobj_tri = copy(testCase.swobj_tri);
            swobj_tri.genmagstr('mode', 'helical', 'S', [1; 0; 0], ...
                                'n', [0 0 1], 'k', [1/3 1/3 0]);
            swobj_tri.genmagstr('mode', 'helical', 'S', [1; 0; 0], ...
                                'n', [0 0 1], 'k', [1/3 1/3 0]);
            expected_mag_str = testCase.default_mag_str;
            expected_mag_str.k = [1/3; 1/3; 0];
            expected_mag_str.F = [1.5; 1.5j; 0];
            testCase.verify_obj(expected_mag_str, swobj_tri.mag_str);

        end
        function test_helical_tri_existing_k(testCase)
            swobj_tri = copy(testCase.swobj_tri);
            swobj_tri.genmagstr('mode', 'helical', 'S', [1; 0; 0], ...
                                'n', [0 0 1], 'k', [1/3 1/3 0]);
            expected_mag_str = testCase.default_mag_str;
            expected_mag_str.k = [1/3; 1/3; 0];
            expected_mag_str.F = [1.5; 1.5j; 0];
            testCase.verify_obj(expected_mag_str, swobj_tri.mag_str);

        end
        function test_direct_afm_chain(testCase)
            afm_chain = copy(testCase.swobj);
            afm_chain.genmagstr('mode', 'direct', 'k',[1/2 0 0], ...
                                'n',[1 0 0],'S',[0 0; 1 -1;0 0], ...
                                'nExt',[2 1 1]);
            expected_mag_str = testCase.default_mag_str;
            expected_mag_str.nExt(1) = int32(2);
            expected_mag_str.k = [0.5; 0; 0];
            expected_mag_str.F = [0 0; 1 -1; 0 0];
            testCase.verify_obj(expected_mag_str, afm_chain.mag_str);

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
    end

end