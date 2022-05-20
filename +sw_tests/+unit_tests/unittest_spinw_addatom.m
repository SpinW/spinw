classdef unittest_spinw_addatom < matlab.mock.TestCase

    properties
        swobj = [];
    end

    methods (TestMethodSetup)
        function setup_spinw_model(testCase)
            testCase.swobj = spinw(); % default init
        end
    end

    methods (Test)
        function test_no_input_calls_help(testCase)
            % Tests that if call addmatrix with no input, it calls the help
            help_function = sw_tests.utilities.mock_function('swhelp');
            testCase.swobj.addatom();
            testCase.assertEqual(help_function.n_calls, 1);
            testCase.assertEqual(help_function.arguments, ...
                {{'spinw.addatom'}});
        end
        
        function test_add_multiple_atom_throws_mismatch_spin(testCase)
            testCase.verifyError(...
                @() testCase.swobj.addatom('r', [0 0;0 0.5;0 0.5], 'S', 1), ...
                'spinw:sw_valid:SizeMismatch')
        end
        
        function test_add_atom_fails_invalid_position_size(testCase)
            testCase.verifyError(...
                @() testCase.swobj.addatom('r', [0, 0]), ...
                'MATLAB:NonIntegerInput')  % need a better error for this
        end
        
        function test_add_atom_fails_without_position(testCase)
            testCase.verifyError(...
                @() testCase.swobj.addatom('S', 1), ...
                'sw_readparam:MissingParameter')
        end
        
        function test_add_atom_fails_with_negative_spin(testCase)
            testCase.verifyError(...
                @() testCase.swobj.addatom('r', [0; 0; 0], 'S', -1), ...
                'spinw:addatom:WrongInput')
        end
        
        function test_add_single_default_atom_with_only_position(testCase)
            pos = [0; 0; 0];
            testCase.swobj.addatom('r', pos)
            unit_cell = testCase.swobj.unit_cell;
            testCase.assertEqual(unit_cell.r, pos)
            testCase.assertEqual(unit_cell.S, 0)  % default non-mag
            testCase.assertEqual(unit_cell.A, int32(-1))  % natural isotope mix
            testCase.assertEqual(unit_cell.Z, int32(113))  % same for all subsequent default atoms added
            testCase.assertEqual(unit_cell.ox, 0)
            testCase.assertEqual(unit_cell.occ, 1)
            testCase.assertEqual(unit_cell.label, {'atom_1'})
            testCase.assertEqual(unit_cell.color, int32([255; 0; 0]))  % same for all subsequent default atoms added
        end
        
        function test_add_multiple_atom_with_only_position(testCase)
            pos = [[0; 0; 0] [0; 0; 0.5]];
            testCase.swobj.addatom('r', pos)
            unit_cell = testCase.swobj.unit_cell;
            testCase.assertEqual(unit_cell.r, pos)
            testCase.assertEqual(unit_cell.S, [0,0])  % default non-mag
            testCase.assertEqual(unit_cell.label, {'atom_1', 'atom_2'})
        end

        
        function test_add_atom_with_update_true_different_spin(testCase)
            pos = [0; 0; 0];
            testCase.swobj.addatom('r', pos, 'label', 'atom1')
            testCase.swobj.addatom('r', pos, 'S', 1, 'label', 'atom1')
            unit_cell = testCase.swobj.unit_cell;
            testCase.assertEqual(size(unit_cell.r), [3, 1])  % 1 atom
            testCase.assertEqual(unit_cell.S, 1)  % updated spin
        end
        
        function test_add_atom_update_false_different_spin(testCase)
            pos = [0; 0; 0];
            testCase.swobj.addatom('r', pos, 'label', 'atom1')
            testCase.swobj.addatom('r', pos, 'S', 1, 'label', 'atom1', 'update', false)
            unit_cell = testCase.swobj.unit_cell;
            testCase.assertEqual(size(unit_cell.r), [3, 2])  % 2 atoms
            testCase.assertEqual(unit_cell.S, [0, 1])
        end
        
        function test_add_atom_will_not_update_different_pos(testCase)
            pos = [0; 0; 0];
            testCase.swobj.addatom('r', pos, 'label', 'atom1')
            testCase.swobj.addatom('r', pos + 0.5, 'S', 1, 'label', 'atom1')
            unit_cell = testCase.swobj.unit_cell;
            testCase.assertEqual(size(unit_cell.r), [3, 2])  % 2 atom
            testCase.assertEqual(unit_cell.S, [0, 1])
        end
        
        function test_add_atom_named_ion_lookup(testCase)
            testCase.swobj.addatom('r', [0; 0; 0], 'label', 'Mn3+')
            unit_cell = testCase.swobj.unit_cell;
            testCase.assertEqual(size(unit_cell.r), [3, 1])  % 2 atom
            testCase.assertEqual(unit_cell.S, 2)
            testCase.assertEqual(unit_cell.ox, 3)
            testCase.assertEqual(unit_cell.Z, int32(25))
            testCase.assertEqual(unit_cell.b(1), -3.73) % coh scat. len
        end
        
     end

end
