classdef unittest_spinw_addcoupling < sw_tests.unit_tests.unittest_super

    properties
        swobj = [];
        default_coupling = struct('dl', int32([1, 0, 0, 1, 1, 1, 0, 1, 0;
                                               0, 1, 0, 0,-1, 1,-1, 0, 1;
                                               0, 0, 1,-1, 0, 0, 1, 1, 1]), ...
            'atom1', ones(1,9, 'int32'), 'atom2', ones(1, 9, 'int32'), ...
            'mat_idx', zeros(3, 9, 'int32'), 'idx', int32([1 1 1 2 2 2 2 2 2]), ...
            'type', zeros(3, 9, 'int32'), 'sym', zeros(3, 9, 'int32'), ...
            'rdip', 0.0, 'nsym',int32(0))
    end
    properties (TestParameter)
        bond_atoms = {{'atom_1', 'atom_1'}, 'atom_1', [1,1], 1};
        invalid_bond_atoms = {{'atom_1', 'atom_2', 'atom_1'}, [1,2,1]}
        mat_label = {'J1', 1}
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
            testCase.verifyError(...
                @() testCase.swobj.addcoupling('mat', 'J2', 'bond', 1), ...
                'spinw:addcoupling:WrongMatrixLabel')
        end
        
        function test_add_coupling_to_sym_equivalent_bonds(testCase, mat_label)
            testCase.swobj.addcoupling('mat', mat_label, 'bond', 1)
            % check matrix added to first three (symm equiv.) bonds
            expected_coupling = testCase.default_coupling;
            expected_coupling.mat_idx(1,1:3) = 1;
            expected_coupling.sym(1,1:3) = 1;
            testCase.verify_val(expected_coupling, ...
                testCase.swobj.coupling)
        end
        
        function test_add_coupling_to_individual_bond(testCase)
            testCase.swobj.addcoupling('mat', 'J1', 'bond', 1, 'subIdx', 1)
            % check matrix added to only first bond with subIdx = 1
            expected_coupling = testCase.default_coupling;
            expected_coupling.mat_idx(1,1) = 1;
            testCase.verify_val(expected_coupling, ...
                testCase.swobj.coupling)
        end
        
        function test_add_ccoupling_lower_symm_with_subIdx(testCase)
            testCase.swobj.genlattice('spgr', 'I 4')
            testCase.swobj.gencoupling('maxDistance',5) % generate bond list
            testCase.verifyWarning(...
                @() testCase.swobj.addcoupling('mat', 'J1', 'bond', 1, 'subIdx', 1), ...
                'spinw:addcoupling:SymetryLowered')
            testCase.assertFalse(any(testCase.swobj.coupling.sym(1,:)))
        end
        
        function test_add_coupling_uses_subIdx_only_first_bond(testCase)
            testCase.verifyWarning(...
                @() testCase.swobj.addcoupling('mat', 'J1', 'bond', [1, 2], 'subIdx', 1), ...
                'spinw:addcoupling:CouplingSize')
            expected_coupling = testCase.default_coupling;
            expected_coupling.mat_idx(1,1) = 1;
            testCase.verify_val(expected_coupling, ...
                testCase.swobj.coupling)
        end
        
        function test_add_coupling_to_multiple_bonds(testCase)
            testCase.swobj.addcoupling('mat', 'J1', 'bond', [1, 2])
            expected_coupling = testCase.default_coupling;
            expected_coupling.mat_idx(1,:) = 1;
            expected_coupling.sym(1,:) = 1;
            testCase.verify_val(expected_coupling, ...
                testCase.swobj.coupling)
        end
        
        function test_add_coupling_to_bonds_using_atom_label(testCase, bond_atoms)
            % add other atom to the sw object made on setup
            testCase.swobj.addatom('r',[0.5, 0.5, 0.5],...
                'S',1, 'label', {'atom_2'})
            testCase.swobj.gencoupling('maxDistance',3) % generate bond list
            % add bond only between atom_1 and atom_1 (bond 2 in
            % obj.coupling.idx)
            testCase.swobj.addcoupling('mat', 'J1', 'bond', 2, ...
                'atom', bond_atoms)
            % check matrix added to atom_1-atom_1 bonds not atom_1-atom_2    
            coupl = testCase.swobj.coupling;
            ibond = find(coupl.mat_idx(1,:));
            testCase.assertEqual(ibond, 9:11);
            testCase.assertTrue(all(coupl.atom1(ibond)==int32(1)))
            testCase.assertTrue(all(coupl.atom2(ibond)==int32(1)))
            testCase.assertTrue(all(coupl.idx(ibond)==int32(2)))
        end
        
        function test_add_coupling_to_bond_between_atoms_different_to_label(testCase, bond_atoms)
            % add other atom to the sw object made on setup
            testCase.swobj.addatom('r',[0.5, 0.5, 0.5],...
                'S',1, 'label', {'atom_2'})
            testCase.swobj.gencoupling('maxDistance',3) % generate bond list
            % atom_1-atom_1 bonds have 'bond' index 2 (in sw.coupling.idx)
            testCase.verifyError(...
                @() testCase.swobj.addcoupling('mat', 'J1', 'bond', 1, ...
                    'atom', bond_atoms), ...
                'spinw:addcoupling:NoBond') 
        end
        
        function test_add_coupling_to_bond_with_invalid_atom_label(testCase)
            testCase.verifyError(...
                @() testCase.swobj.addcoupling('mat', 'J1', 'bond', 1, ...
                'atom', {'atom_1', 'atom_3'}), ...
                'spinw:addcoupling:WrongInput')
        end
        
        function test_add_coupling_more_than_two_atom_labels(testCase, invalid_bond_atoms)
            testCase.verifyError(...
                @() testCase.swobj.addcoupling('mat', 'J1', 'bond', 1, ...
                'atom', invalid_bond_atoms), ...
                'spinw:addcoupling:WrongInput')
        end
        
        function test_add_coupling_biquadratic_exchange(testCase)
            testCase.swobj.addcoupling('mat', 'J1', 'bond', 1, ...
                'type', 'biquadratic')
            expected_coupling = testCase.default_coupling;
            expected_coupling.mat_idx(1,1:3) = 1;
            expected_coupling.sym(1,1:3) = 1;
            expected_coupling.type(1,1:3) = 1;
            testCase.verify_val(expected_coupling, ...
                testCase.swobj.coupling)
        end
        
        function test_add_coupling_invalid_exchange_type(testCase)
            testCase.verifyError(...
                @() testCase.swobj.addcoupling('mat', 'J1', 'bond', 1, ...
                    'type', 'invalid_type'), ...
                'spinw:addcoupling:WrongInput')
        end
        
        function test_add_coupling_multiple_on_same_bond(testCase)
            testCase.swobj.addmatrix('label','J2','value', 1)
            testCase.swobj.addcoupling('mat', 'J1', 'bond', 1)
            testCase.swobj.addcoupling('mat', 2, 'bond', 1)
            % check matrix added to first three (symm equiv.) bonds 
            expected_coupling = testCase.default_coupling;
            expected_coupling.mat_idx(1,1:3) = 1;
            expected_coupling.mat_idx(2,1:3) = 2;
            expected_coupling.sym(1:2,1:3) = 1;
            testCase.verify_val(expected_coupling, ...
                testCase.swobj.coupling)
        end
        
        function test_add_coupling_max_num_matrices_added(testCase)
            % add 4 matrices to a bond
            for imat = 1:4
                mat_str = ['J', num2str(imat)];
                testCase.swobj.addmatrix('label', mat_str,'value', imat)
                if imat < 4
                    testCase.swobj.addcoupling('mat', mat_str, 'bond', 1)
                else
                    % exceeded max. of 3 matrices on a bond
                    testCase.verifyError(...
                        @() testCase.swobj.addcoupling('mat', mat_str, 'bond', 1), ...
                        'spinw:addcoupling:TooManyCoupling')
                    
                end
            end
        end
        
        function test_add_coupling_with_sym_false(testCase)
            testCase.swobj.addcoupling('mat', 'J1', 'bond', 1, 'sym', false)
            expected_coupling = testCase.default_coupling;
            expected_coupling.mat_idx(1,1:3) = 1;
            testCase.verify_val(expected_coupling, ...
                testCase.swobj.coupling)
        end
        
        function test_add_coupling_duplicate_matrix_on_same_bond(testCase)
            testCase.swobj.addcoupling('mat', 'J1', 'bond', 1)
            testCase.verifyWarning(...
                @() testCase.swobj.addcoupling('mat', 'J1', 'bond', 1), ...
                'spinw:addcoupling:CouplingIdxWarning')
            % check matrix only added once
            expected_coupling = testCase.default_coupling;
            expected_coupling.mat_idx(1,1:3) = 1;
            expected_coupling.sym(1,1:3) = 1;
            testCase.verify_val(expected_coupling, ...
                testCase.swobj.coupling)
        end
        
     end

end
