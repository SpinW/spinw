classdef unittest_spinw_genmagstr < sw_tests.unit_tests.unittest_super

    properties
        swobj = [];
        swobj_tri = [];
        default_magstr = struct('S', [0; 1; 0], ...
                                'k', [0 0 0], ...
                                'n', [0 0 1], ...
                                'N_ext', [1 1 1], ...
                                'exact', true);
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
            swobj.genmagstr('mode','direct', 'k',[0 0 0],'n',[1 0 0],'S',[0; 1; 0]);
            testCase.verify_obj(testCase.default_magstr, swobj.magstr);

        end
        function test_direct_fm_chain_nok(testCase)
            swobj = copy(testCase.swobj);
            swobj.genmagstr('mode','direct', 'n',[1 0 0],'S',[0; 1; 0]);
            testCase.verify_obj(testCase.default_magstr, swobj.magstr);
        end
        function test_helical_tri(testCase)
            swobj_tri = copy(testCase.swobj_tri);
            swobj_tri.genmagstr('mode', 'helical', 'S', [1; 0; 0], ...
                                         'n', [0 0 1], 'k', [1/3 1/3 0]);
            expected_magstr = testCase.default_magstr;
            expected_magstr.S = [1.5; 0; 0];
            expected_magstr.k = [1/3 1/3 0];
            testCase.verify_obj(expected_magstr, swobj_tri.magstr);

        end
        function test_direct_afm_chain(testCase)
            afm_chain = copy(testCase.swobj);
            afm_chain.genmagstr('mode', 'direct', 'k',[1/2 0 0], ...
                                'n',[1 0 0],'S',[0 0; 1 -1;0 0], ...
                                'nExt',[2 1 1]);
            expected_magstr = testCase.default_magstr;
            expected_magstr.S = [0 0; 1 -1; 0 0];
            expected_magstr.N_ext = [2 1 1];
            testCase.verify_obj(expected_magstr, afm_chain.magstr);

        end
        function test_random_structure(testCase)
            swobj = copy(testCase.swobj);
            swobj.genmagstr('mode','random');
            magstr1 = swobj.magstr;
            swobj.genmagstr('mode','random');
            magstr2 = swobj.magstr;
            % Check structure is random each time - S is different
            testCase.verifyNotEqual(magstr1.S, magstr2.S);
            % Check S has correct magnitude
            testCase.verify_val(norm(magstr1.S, 2), 1);
            % Check other fields
            expected_magstr = testCase.default_magstr;
            testCase.verify_obj(rmfield(expected_magstr, 'S'), ...
                                rmfield(magstr1, 'S'));
        end
        function test_random_structure_multiatom_and_nExt(testCase)
            swobj = copy(testCase.swobj);
            swobj.addatom('r', [0.5 0.5 0],'S', 2);
            nExt = [2 2 2];
            swobj.genmagstr('mode','random', 'nExt', nExt);
            magstr1 = swobj.magstr;
            swobj.genmagstr('mode','random', 'nExt', nExt);
            magstr2 = swobj.magstr;
            % Check structure is random each time - S is different
            testCase.verifyNotEqual(magstr1.S, magstr2.S);
            % Check S has correct magnitudes
            testCase.verify_val(vecnorm(magstr1.S(:, 1:2:end), 2), ones(1, 8));
            testCase.verify_val(vecnorm(magstr1.S(:, 2:2:end), 2), 2*ones(1, 8));
            % Check other fields
            expected_magstr = testCase.default_magstr;
            expected_magstr.N_ext = nExt;
            testCase.verify_obj(rmfield(expected_magstr, 'S'), ...
                                rmfield(magstr1, 'S'));
        end
    end

end