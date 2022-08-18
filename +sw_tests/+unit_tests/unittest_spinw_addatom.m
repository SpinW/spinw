classdef unittest_spinw_addatom < sw_tests.unit_tests.unittest_super

    properties
        swobj = [];
        default_unit_cell = struct('r', [0; 0; 0], 'S', 0, ...
            'label', {{'atom_1'}}, 'color', int32([255; 0; 0]), 'ox', 0, ...
            'occ', 1, 'b', [1; 1], 'ff', [zeros(2,10) [1; 1]], ...
            'A', int32(-1), 'Z', int32(113), 'biso', 0)
        ff = [0.4198, 14.2829, 0.6054, 5.4689, 0.9241, -0.0088, ...
                0,      0,        0,      0,      -0.9498;
              6.9270, 0.3783, 2.0813, 0.0151, 11.1284, 5.3800, ...
                2.3751, 14.4296, -0.4193, 0.0049, -0.0937]; % Mn3+ form-fac
    end
    properties (TestParameter)
        property_error = {{'Z', 'spinw:addatom:WrongInput'}, ...
            {'A', 'spinw:addatom:WrongInput'}, ...
            {'biso', 'spinw:addatom:WrongInput'}, ...
            {'ox', 'spinw:addatom:WrongInput'}, ...
            {'S', 'spinw:sw_valid:SizeMismatch'}};
        oxidation_label = {{3, 'Fe3+_1'}, {-3, 'Fe3-_1'}};
        pos_vector = {[0;0;0], [0 0 0]}
        property_value = {{'S',1}, {'occ', 0.5}, {'biso', 0.5}};
        b_name = {'b', 'bn'}
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
               
        function test_add_multiple_atom_throws_mismatch_parameter_size(testCase, property_error)
            [prop, error] = property_error{:};
            testCase.verifyError(...
                @() testCase.swobj.addatom('r', [0 0;0 0.5;0 0.5], prop, 2), ...
                error)
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
        
        function test_add_atom_warns_bx_provided(testCase)
            testCase.verifyWarning(...
                @() testCase.swobj.addatom('r', [0; 0; 0], 'bx', 2), ...
                'spinw:addatom:DeprecationWarning')
        end
        
        function test_add_atom_with_bn_throws_deprecation_warning(testCase)
            testCase.verifyWarning(...
                @() testCase.swobj.addatom('r', [0; 0; 0], 'bn', 2), ...
                'spinw:addatom:DeprecationWarning')
        end
        
        function test_add_atom_warns_b_and_bn_provided(testCase)
            testCase.verifyWarning(...
                @() testCase.swobj.addatom('r', [0; 0; 0], 'b', 2, 'bn', 3), ...
                'spinw:addatom:WrongInput')
            % check value for scattering length is bn provided
            testCase.assertEqual(testCase.swobj.unit_cell.b(1,1), 3)
        end
        
        function test_add_single_default_atom_with_only_position(testCase, pos_vector)
            testCase.swobj.addatom('r', pos_vector)
            testCase.verify_val(testCase.swobj.unit_cell, ...
                testCase.default_unit_cell);
        end
        
        function test_add_single_atom_custom_parameters(testCase, property_value)
            [prop, val] = property_value{:};
            testCase.swobj.addatom('r', [0; 0; 0], prop, val)
            expected_unit_cell = testCase.default_unit_cell;
            expected_unit_cell.(prop) = val;
            testCase.verify_val(testCase.swobj.unit_cell, ...
                expected_unit_cell)
            
        end
        
        function test_add_atom_with_custom_scatt_length(testCase, b_name)
            b = 2;
            testCase.swobj.addatom('r', [0;0;0], b_name, b)
            expected_unit_cell = testCase.default_unit_cell;
            expected_unit_cell.b = [b; 1];
            testCase.verify_val(testCase.swobj.unit_cell, ...
                expected_unit_cell)
        end
        
        function test_add_multiple_atom_with_single_call(testCase)
            pos = [[0; 0; 0] [0; 0; 0.5]];
            S = [0, 1];
            testCase.swobj.addatom('r', pos, 'S', S)
            unit_cell = testCase.swobj.unit_cell;
            testCase.assertEqual(unit_cell.r, pos)
            testCase.assertEqual(unit_cell.S, S)  % default non-mag
            testCase.assertEqual(unit_cell.label, {'atom_1', 'atom_2'})
        end

        
        function test_add_atom_with_update_true_different_spin(testCase)
            pos = [0; 0; 0];
            label = 'atom';
            testCase.swobj.addatom('r', pos, 'label', label)
            testCase.swobj.addatom('r', pos, 'S', 1, 'label', label)
            expected_unit_cell = testCase.default_unit_cell;
            expected_unit_cell.S = 1;
            expected_unit_cell.label = {label};
            testCase.verify_val(testCase.swobj.unit_cell, ...
                expected_unit_cell)
        end
        
        function test_add_atom_update_false_different_spin(testCase)
            pos = [0; 0; 0];
            testCase.swobj.addatom('r', pos, 'label', 'atom1')
            testCase.verifyWarning(...
                @() testCase.swobj.addatom('r', pos, 'S', 1, ...
                    'label', 'atom1', 'update', false), ...
                'spinw:addatom:WrongInput')  % warns occ > 1
            unit_cell = testCase.swobj.unit_cell;
            testCase.assertEqual(unit_cell.r, [pos, pos])  % 2 atoms
            testCase.assertEqual(unit_cell.S, [0, 1])
        end
        
        function test_add_atom_will_not_update_different_pos(testCase)
            pos = [0; 0; 0];
            testCase.swobj.addatom('r', pos, 'label', 'atom1')
            testCase.swobj.addatom('r', pos + 0.5, 'S', 1, 'label', 'atom1')
            unit_cell = testCase.swobj.unit_cell;
            testCase.assertEqual(unit_cell.r, [pos, pos+0.5])  % 2 atom
            testCase.assertEqual(unit_cell.S, [0, 1])
        end
        
        function test_add_atom_named_ion_lookup(testCase)
            label = 'Mn3+';
            testCase.swobj.addatom('r', [0; 0; 0], 'label', label)
            expected_unit_cell = testCase.default_unit_cell;
            expected_unit_cell.S = 2;
            expected_unit_cell.ox = 3;
            expected_unit_cell.Z = int32(25);
            expected_unit_cell.b = [-3.73; 1];
            expected_unit_cell.ff = testCase.ff;
            expected_unit_cell.label = {label};
            expected_unit_cell.color = int32([156; 122; 199]);
            testCase.verify_val(testCase.swobj.unit_cell, ...
                expected_unit_cell, 'abs_tol', 1e-4)
        end
        
        function test_add_atom_lookup_by_Z_with_custom_ox(testCase, oxidation_label)
            [ox, label] = oxidation_label{:};
            Z = int32(26);
            testCase.swobj.addatom('r',[0;0;0], 'Z', Z, 'ox', ox)
            expected_unit_cell = testCase.default_unit_cell;
            expected_unit_cell.ox = ox;
            expected_unit_cell.Z = Z;
            expected_unit_cell.label = {label};
            expected_unit_cell.color = int32([224; 102; 51]);
            testCase.verify_val(testCase.swobj.unit_cell, ...
                expected_unit_cell)
        end
        
        function test_add_atom_with_custom_form_factor(testCase)
            label = 'Mn3+';
            testCase.swobj.addatom('r', [0; 0; 0], 'label', label, ...
                'formfact', 1:9);
            expected_unit_cell = testCase.default_unit_cell;
            expected_unit_cell.ox = 3;
            expected_unit_cell.Z = int32(25);
            expected_unit_cell.b = [-3.73; 1];
            expected_unit_cell.ff = [1:8 0 0 9; testCase.ff(2,:)];
            expected_unit_cell.label = {label};
            expected_unit_cell.color = int32([156; 122; 199]);
            testCase.verify_val(testCase.swobj.unit_cell, ...
                expected_unit_cell, 'abs_tol', 1e-4)
        end
        
        function test_add_atom_with_custom_form_factor_wrong_size(testCase)
            testCase.verifyError(...
                @() testCase.swobj.addatom('r', [0; 0; 0], 'label', 'Mn3+', ...
                    'formfact', 1:8), ...
                'MATLAB:catenate:dimensionMismatch');
        end    
                
    end
     
    methods (Test, TestTags = {'Symbolic'})
        function test_add_atom_in_symbolic_mode_has_symbolic_S(testCase)
            pos = [0 0.5; 0 0.5; 0 0.5];
            testCase.swobj.symbolic(true)
            testCase.swobj.addatom('r', pos, 'S', [0,1])
            unit_cell = testCase.swobj.unit_cell;
            testCase.assertEqual(unit_cell.r, pos)
            testCase.assertEqual(unit_cell.S, [sym(0), sym('S_2')])
        end
    end

end
