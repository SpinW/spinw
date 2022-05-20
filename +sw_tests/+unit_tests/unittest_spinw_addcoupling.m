classdef unittest_spinw_addcoupling < matlab.mock.TestCase

    properties
        swobj = [];
    end

    methods (TestMethodSetup)
        function setup_spinw_model(testCase)
            testCase.swobj = spinw(); % default init
            testCase.swobj.addatom('r',[0 0 0],'S',1)
            testCase.swobj.gencoupling('maxDistance',5) % generate bond list
            testCase.swobj.addmatrix('label','J1','value', 1)
        end
    end

    methods (Test)
            
        function test_add_coupling_requires_mat_and_bond(testCase)
            testCase.verifyError(...
                @() testCase.swobj.addcoupling(), ...
                'sw_readparam:MissingParameter')
            testCase.verifyError(...
                @() testCase.swobj.addcoupling('mat','J1'), ...
                'sw_readparam:MissingParameter')
            testCase.verifyError(...
                @() testCase.swobj.addcoupling('bond', 1), ...
                'sw_readparam:MissingParameter')
        end
        
        function test_add_coupling_requires_run_gencoupling(testCase)
            sw = spinw(); % default init
            sw.addatom('r',[0 0 0],'S',1)
            sw.addmatrix('label','J1','value', 1)
            testCase.verifyError(...
                @() sw.addcoupling('mat','J1','bond',1), ...
                'spinw:addcoupling:CouplingError')
        end
        
        function test_add_coupling_requires_magnetic_atom(testCase)
            sw = spinw(); % default init
            sw.addatom('r',[0 0 0],'S',0)
            testCase.verifyError(...
                @() sw.addcoupling('mat','J1','bond',1), ...
                'spinw:addcoupling:NoMagAtom')
        end
        
        function test_add_coupling_with_invalid_matrix(testCase)
            testCase.verifyWarning(...
                @() testCase.swobj.addcoupling('mat', 'J2', 'bond', 1), ...
                'spinw:addcoupling:WrongMatrixLabel')
        end
        
        function test_add_coupling_to_sym_equivalent_bonds(testCase)
            testCase.swobj.addcoupling('mat', 'J1', 'bond', 1)
            % check matrix added to first three (symm equiv.) bonds 
            coupl = testCase.swobj.coupling;
            testCase.assertEqual(coupl.mat_idx(1,1:3), int32([1 ,1, 1]))
            testCase.assertFalse(any(coupl.mat_idx(1,4:end)))
            testCase.assertTrue(all(coupl.sym(1,1:3)))
            % check default quadratic exchange
            testCase.assertFalse(any(testCase.swobj.coupling.type(1,:)))
        end
        
        function test_add_coupling_using_matrix_index(testCase)
            testCase.swobj.addcoupling('mat', 1, 'bond', 1)
            testCase.assertEqual(testCase.swobj.coupling.mat_idx(1,1:3), ...
                int32([1 ,1, 1]))
        end
        
        function test_add_coupling_to_individual_bond(testCase)
            testCase.swobj.addcoupling('mat', 'J1', 'bond', 1, 'subIdx', 1)
            % check matrix added to only first bond with subIdx = 1
            coupl = testCase.swobj.coupling;
            testCase.assertEqual(coupl.mat_idx(1,1), int32(1))
            testCase.assertFalse(any(coupl.mat_idx(1,2:end)))
            % check removes symmetry of that bond
            testCase.assertFalse(any(coupl.sym(1,1)))
        end
        
        function test_add_coupling_to_multiple_bonds(testCase)
            testCase.swobj.addcoupling('mat', 'J1', 'bond', [1, 2])
            coupl = testCase.swobj.coupling;
            testCase.assertTrue(all(coupl.mat_idx(1,:)==int32(1)))
        end
        
        function test_add_coupling_to_bonds_using__atom_label(testCase)
            % add other atom to the sw object made on setup
            testCase.swobj.addatom('r',[0.5, 0; 0 0.5; 0 0],...
                'S',[1,1], 'label', {'atom_1', 'atom_2'})
            testCase.swobj.gencoupling('maxDistance',2) % generate bond list
            % add bond only between atom_1 and atom_2
            testCase.swobj.addcoupling('mat', 'J1', 'bond', 1, ...
                'atom', {'atom_1', 'atom_2'})
            % check matrix added to atom_1-atom_2 bonds not atom_1-atom_1
            testCase.assertEqual(testCase.swobj.coupling.mat_idx(1,:), ...
                int32([0, 1, 0, 1]))
        end
        
        function test_add_coupling_biquadratic_exchange(testCase)
            testCase.swobj.addcoupling('mat', 'J1', 'bond', 1, 'type', 'biquadratic')
            testCase.assertEqual(testCase.swobj.coupling.type(1,1:3), ...
                int32([1 ,1, 1]))
        end
        
        function test_add_coupling_multiple_on_same_bond(testCase)
            testCase.swobj.addmatrix('label','J2','value', 1)
            testCase.swobj.addcoupling('mat', 'J1', 'bond', 1)
            testCase.swobj.addcoupling('mat', 2, 'bond', 1)
            % check matrix added to first three (symm equiv.) bonds 
            coupl = testCase.swobj.coupling;
            testCase.assertEqual(coupl.mat_idx(1,1:3), int32([1, 1, 1]))
            testCase.assertEqual(coupl.mat_idx(2,1:3), int32([2, 2, 2]))
        end
        
     end

end
