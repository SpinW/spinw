classdef unittest_spinw_gencoupling < matlab.mock.TestCase

    properties
        swobj = [];
    end

    methods (TestMethodSetup)
        function setup_spinw_model(testCase)
            testCase.swobj = spinw(); % default init
            testCase.swobj.genlattice('lat_const',[3 3 5])
            testCase.swobj.addatom('r',[0 0 0],'S',1)
            testCase.swobj.gencoupling('maxDistance',5) % generate bond list
        end
    end

    methods (Test)
            
        function test_gencoupling_requires_magnetic_atom(testCase)
            % change previously added atom to have S=0 (non magnetic)
            testCase.swobj.addatom('r',[0 0 0],'S',0, ...
                'label', 'atom_1', 'update', true)
            testCase.verifyError(...
                @() testCase.swobj.gencoupling(), ...
                'spinw:gencoupling:NoMagAtom')
        end
        
        function test_gencoupling_with_maxDistance_less_than_dMin(testCase)
            testCase.verifyError(...
                @() testCase.swobj.gencoupling('maxDistance', 0.1, 'dMin', 0.5), ...
                'spinw:gencoupling:MaxDLessThanMinD')
        end
        
        function test_gencoupling_with_maxDistance_less_than_lattice_param(testCase)
            testCase.verifyError(...
                @() testCase.swobj.gencoupling('maxDistance', ...
                testCase.swobj.lattice.lat_const(1)/2), ...
                'spinw:gencoupling:maxDistance')
        end
        
        function test_gencoupling_with_tol_when_maxDistance_equal_lattice_param(testCase)
            % check bonds are created for atoms separated by latt. param.
            % when tol > delta
            delta = 1e-2;
            testCase.swobj.gencoupling('tol', 10*delta, 'maxDistance', ...
                testCase.swobj.lattice.lat_const(1) - delta)
            testCase.assertEqual(testCase.swobj.coupling.atom1, int32([1,1]))
            testCase.assertEqual(testCase.swobj.coupling.atom2, int32([1,1]))
            testCase.assertEqual(testCase.swobj.coupling.dl, ...
                int32([1, 0; 0 1; 0 0]))
        end
        
        function test_gencoupling_with_tolDist_for_symm_equiv_bonds(testCase)
            % change lattice param to slightly break tetragonal sym.
            delta = 1e-3;
            testCase.swobj.genlattice('lat_const',[3 3+delta 5])
            % check that when tolDist > delta the bonds are equiv.
            testCase.swobj.gencoupling('maxDistance', 4, 'tolDist', 10*delta)
            testCase.assertEqual(testCase.swobj.coupling.idx, int32([1,1]))
            % check that when tolDist > delta the bonds are inequiv.
            testCase.swobj.gencoupling('maxDistance', 4, 'tolDist', 0.1*delta)
            testCase.assertEqual(testCase.swobj.coupling.idx, int32([1,2]))
        end
        
        function test_gencoupling_with_non_P0_spacegroup(testCase)
            testCase.swobj.genlattice('spgr', 'P 4')  % not overwrite abc
            testCase.swobj.gencoupling('maxDistance', 4)
            testCase.assertEqual(testCase.swobj.coupling.nsym, int32(1))
        end
        
        % also validate maxSym < dMin
        
     end

end
