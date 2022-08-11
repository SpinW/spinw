classdef unittest_spinw_addaniso < sw_tests.unit_tests.unittest_super

    properties
        swobj = [];
        default_single_ion = struct('aniso', int32(1), ...
            'g', zeros(1,0,'int32'), 'field', [0,0,0], 'T', 0)
    end

    methods (TestMethodSetup)
        function setup_spinw_model(testCase)
            testCase.swobj = spinw();
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
        
        function test_addaniso_with_wrong_atom_label(testCase)
            testCase.verifyError(...
                @() testCase.swobj.addaniso('A1', 'atom_2'), ...
                'spinw:addaniso:WrongString')
        end
        
        function test_addaniso_with_no_magnetic_atom(testCase)
            testCase.swobj.addatom('r',[0,0,0], 'S',0, ...
                'label', 'atom_1', 'update', true)
            testCase.verifyError(...
                @() testCase.swobj.addaniso('A1'), ...
                'spinw:addaniso:NoMagAtom')
        end
        
        function test_addaniso_all_symm_equiv_atoms(testCase)
            testCase.swobj.genlattice('sym','I 4'); % body-centred
            testCase.swobj.addaniso('A1')
            expected_single_ion = testCase.default_single_ion;
            expected_single_ion.aniso = int32([1, 1]);
            testCase.verify_val(expected_single_ion, ...
                testCase.swobj.single_ion)
        end
        
        function test_addaniso_specific_atoms_wrong_atomIdx(testCase)
            testCase.verifyError(...
                @() testCase.swobj.addaniso('A1', 'atom_1', 3), ...
                'MATLAB:matrix:singleSubscriptNumelMismatch')
        end
        
        function test_addaniso_with_atomIdx_error_when_high_symm(testCase)
            testCase.swobj.genlattice('sym','I 4'); % body-centred
            testCase.verifyError(...
                @() testCase.swobj.addaniso('A1', 'atom_1', 1), ...
                'spinw:addaniso:SymmetryProblem')
        end
        
        function test_addaniso_specific_atom_label(testCase)
            % setup unit cell with body-centred atom as different species
            testCase.swobj.addatom('r',[0.5; 0.5; 0.5], 'S',1)
            testCase.swobj.addaniso('A1', 'atom_2')
            expected_single_ion = testCase.default_single_ion;
            expected_single_ion.aniso = int32([0, 1]);
            testCase.verify_val(expected_single_ion, ...
                testCase.swobj.single_ion)
        end
        
        function test_addaniso_overwrites_previous_aniso(testCase)
            testCase.swobj.addmatrix('label','A2','value',diag([-0.5 0 0]))
            testCase.swobj.addaniso('A1')
            testCase.swobj.addaniso('A2')
            expected_single_ion = testCase.default_single_ion;
            expected_single_ion.aniso = int32(2);
            testCase.verify_val(expected_single_ion, ...
                testCase.swobj.single_ion)
        end
        
     end

end
