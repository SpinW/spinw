classdef unittest_spinw_addg < sw_tests.unit_tests.unittest_super

    properties
        swobj = [];
        default_single_ion = struct('g', int32(1), ...
            'aniso', zeros(1,0,'int32'), 'field', [0,0,0], 'T', 0)
    end

    methods (TestMethodSetup)
        function setup_spinw_model(testCase)
            testCase.swobj = spinw();
            testCase.swobj.addatom('r',[0,0,0], 'S',1)
            testCase.swobj.addmatrix('label','g1','value',diag([2 1 1]))
        end
    end

    methods (Test)
            
        function test_addg_requires_matrix(testCase)
            testCase.verifyError(@() testCase.swobj.addg(), ...
                'MATLAB:minrhs') % better if sw_readparam:MissingParameter
        end
        
        function test_addg_wrong_matrix_label(testCase)
            testCase.verifyError(...
                @() testCase.swobj.addg('g2'), ...
                'spinw:addg:WrongCouplingTypeIdx')
        end
        
        function test_addaniso_with_no_magnetic_atom(testCase)
            testCase.swobj.addatom('r',[0,0,0], 'S',0, ...
                'label', 'atom_1', 'update', true)
            testCase.verifyError(...
                @() testCase.swobj.addg('g1'), ...
                'spinw:addg:NoMagAtom')
        end
        
        function test_addg_validates_gmatrix(testCase)
            testCase.swobj.addmatrix('label', 'g', ...
                'value', ones(3));
            testCase.verifyError(@() testCase.swobj.addg('g'), ...
                'spinw:addg:InvalidgTensor')
        end
        
        function test_addg_with_wrong_atom_label_not_write_g(testCase)
            testCase.swobj.addg('g1', 'atom_2')
            testCase.assertFalse(any(testCase.swobj.single_ion.g))
        end
        
        function test_addg_all_symm_equiv_atoms(testCase)
            testCase.swobj.genlattice('sym','I 4'); % body-centred
            testCase.swobj.addg('g1')
            expected_single_ion = testCase.default_single_ion;
            expected_single_ion.g = int32([1, 1]);
            testCase.verify_val(expected_single_ion, ...
                testCase.swobj.single_ion)
        end
        
        function test_addg_with_atomIdx_error_when_high_symm(testCase)
            testCase.swobj.genlattice('sym','I 4'); % body-centred
            testCase.verifyError(...
                @() testCase.swobj.addg('g1', 'atom_1', 1), ...
                'spinw:addg:SymmetryProblem')
        end
        
        function test_addg_specific_atoms_wrong_atomIdx(testCase)
            testCase.verifyError(...
                @() testCase.swobj.addg('g1', 'atom_1', 3), ...
                'MATLAB:matrix:singleSubscriptNumelMismatch')
        end
        
        function test_addg_specific_atom_label(testCase)
            % setup unit cell with body-centred atom as different species
            testCase.swobj.addatom('r',[0.5; 0.5; 0.5], 'S',1)
            testCase.swobj.addg('g1', 'atom_2')
            expected_single_ion = testCase.default_single_ion;
            expected_single_ion.g = int32([0, 1]);
            testCase.verify_val(expected_single_ion, ...
                testCase.swobj.single_ion)
        end
        
        function test_addg_overwrites_previous_gtensor(testCase)
            testCase.swobj.addmatrix('label','g2','value',diag([1, 1, 2]))
            testCase.swobj.addg('g1')
            testCase.swobj.addg('g2')
            expected_single_ion = testCase.default_single_ion;
            expected_single_ion.g = int32(2);
            testCase.verify_val(expected_single_ion, ...
                testCase.swobj.single_ion)
        end
        
     end

end
