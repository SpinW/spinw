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
            testCase.swobj.gencoupling('maxDistance', 7);
            testCase.swobj.addmatrix('value', -eye(3), 'label', 'Ja');
            testCase.swobj.addcoupling('mat', 'Ja', 'bond', 1);
        end
        function setup_tri_model(testCase)
            % Create a simple triangular lattice model
            testCase.swobj_tri = spinw;
            testCase.swobj_tri.genlattice('lat_const', [4 4 6], 'angled', [90 90 120]);
            testCase.swobj_tri.addatom('r', [0 0 0], 'S', 3/2, 'label', 'MCr3');
            J1 = 1;
            testCase.swobj_tri.addmatrix('label','J1','value',J1);
            testCase.swobj_tri.gencoupling;
            testCase.swobj_tri.addcoupling('mat','J1','bond',1);
        end
    end

    methods (Test)
        function test_no_magnetic_atoms_raises_error(testCase)
            swobj = spinw;
            testCase.verifyError(...
                @() swobj.genmagstr(), ...
                'spinw:genmagstr:NoMagAtom')
        end
        function test_chain(testCase)
            swobj = copy(testCase.swobj);
            swobj.genmagstr('mode','direct', 'k',[0 0 0],'n',[1 0 0],'S',[0; 1; 0]);
            testCase.verify_obj(testCase.default_magstr, swobj.magstr);

        end
        function test_chain_nok(testCase)
            swobj = copy(testCase.swobj);
            swobj.genmagstr('mode','direct', 'n',[1 0 0],'S',[0; 1; 0]);
            testCase.verify_obj(testCase.default_magstr, swobj.magstr);
        end
        function test_tri(testCase)
            swobj_tri = copy(testCase.swobj_tri);
            swobj_tri.genmagstr('mode', 'helical', 'S', [1; 0; 0], ...
                                         'n', [0 0 1], 'k', [1/3 1/3 0]);
            expected_magstr = testCase.default_magstr;
            expected_magstr.S = [1.5; 0; 0];
            expected_magstr.k = [1/3 1/3 0];
            testCase.verify_obj(expected_magstr, swobj_tri.magstr);

        end
    end

end
