classdef unittest_spinw_addaniso < matlab.mock.TestCase

    properties
        swobj = [];
    end

    methods (TestMethodSetup)
        function setup_spinw_model(testCase)
            testCase.swobj = spinw();
            testCase.swobj.genlattice('spgr','I 4'); % body-centred
            testCase.swobj.addatom('r',[0,0,0], 'S',1)
            testCase.swobj.addmatrix('label','A1','value',diag([-0.1 0 0]))
        end
    end

    methods (Test)
            
        function test_addaniso_requires_matrix(testCase)
            testCase.verifyError(...
                @() testCase.swobj.addaniso(), ...
                'MATLAB:minrhs') % better if sw_readparam:MissingParameter
        end
        
        function test_addaniso_wrong_matrix_label(testCase)
            testCase.verifyError(...
                @() testCase.swobj.addaniso('A2'), ...
                'spinw:addaniso:WrongCouplingTypeIdx')
        end
        
        function test_addaniso_with_wrong_atom_label_not_write_g(testCase)
            testCase.swobj.addg('A1', 'atom_2')
            testCase.assertFalse(any(testCase.swobj.single_ion.aniso))
        end
        
        function test_addaniso_all_symm_equiv_atomsg(testCase)
            testCase.swobj.addaniso('A1')
            testCase.assertEqual(testCase.swobj.single_ion.aniso, ...
                int32([1, 1]))
        end
        
        function test_addaniso_specific_atoms_ignoring_symm(testCase)
            testCase.swobj.addaniso('A1', 'atom_1', 1)  % only corner atoms
            testCase.assertEqual(testCase.swobj.single_ion.aniso, int32([1, 0]))
        end
        
        function test_addaniso_specific_atoms_wrong_atom_idx(testCase)
            testCase.verifyError(...
                @() testCase.swobj.addaniso('A1', 'atom_1', 3), ...
                'MATLAB:matrix:singleSubscriptNumelMismatch')
        end
        
        function test_addaniso_specific_atom_label(testCase)
            % add body-centred atom as different species (revert P1 symm)
            testCase.swobj.genlattice('spgr','P 1');
            testCase.swobj.addatom('r',[0.5, 0.5, 0.5], 'S',2)
            testCase.swobj.addaniso('A1', 'atom_2')
            testCase.assertEqual(testCase.swobj.single_ion.aniso, int32([0, 1]))
        end
     end

end
