classdef unittest_spinw_addmatrix < sw_tests.unit_tests.unittest_super

    properties
        swobj = [];
        default_matrix = struct('mat', eye(3), 'color', int32([0;0;0]), ...
            'label', {{'mat1'}})
    end
    properties (TestParameter)
        wrong_matrix_input = {'str', []};
    end

    methods (TestMethodSetup)
        function setup_spinw_model(testCase)
            testCase.swobj = spinw;
        end
    end

    methods (Test)
        function test_no_input_calls_help(testCase)
            % Tests that if call addmatrix with no input, it calls the help
            help_function = sw_tests.utilities.mock_function('swhelp');
            testCase.swobj.addmatrix();
            testCase.assertEqual(help_function.n_calls, 1);
            testCase.assertEqual(help_function.arguments, ...
                {{'spinw.addmatrix'}});
        end
        
        function test_non_numeric_input_warns_and_not_add_matrix(testCase, wrong_matrix_input)
            testCase.verifyError(...
                @() testCase.swobj.addmatrix('value', wrong_matrix_input), ...
                'spinw:addmatrix:WrongInput');
            testCase.assertTrue(isempty(testCase.swobj.matrix.mat));
        end
        
        function test_label_and_no_matrix_warns_and_adds_default(testCase)
            testCase.verifyWarning(...
                @() testCase.swobj.addmatrix('label', 'mat'), ...
                'spinw:addmatrix:NoValue');
            testCase.assertEqual(testCase.swobj.matrix.mat, eye(3));
        end
        
        function test_single_value_adds_diag_matrix_no_label_color(testCase)
            J = 1.0;
            testCase.swobj.addmatrix('value', J);
            expected_matrix = testCase.default_matrix;
            expected_matrix.mat = J*eye(3);
            testCase.verify_spinw_matrix(testCase.swobj.matrix, ...
                                         expected_matrix)
        end
        
        function test_matrix_value_added_no_modification(testCase)
            mat = reshape(1:9, 3, 3);
            testCase.swobj.addmatrix('value', mat);
            expected_matrix = testCase.default_matrix;
            expected_matrix.mat = mat;
            testCase.verify_spinw_matrix(testCase.swobj.matrix, ...
                                         expected_matrix)
        end
        
        function test_vector_value_adds_DM_matrix(testCase)
            [m1, m2, m3] = deal(1,2,3);
            testCase.swobj.addmatrix('value', [m1, m2, m3]);
            expected_matrix = testCase.default_matrix;
            expected_matrix.mat = [0 m3 -m2; -m3 0 m1; m2 -m1 0];
            testCase.verify_spinw_matrix(testCase.swobj.matrix, ...
                                         expected_matrix)
        end
        
        function test_add_multiple_matrices_separate_calls(testCase)
            nmat = 2;
            for imat = 1:nmat
                testCase.swobj.addmatrix('value', imat);
            end
            % check matrix properties have correct dimensions
            expected_matrix = testCase.default_matrix;
            expected_matrix.mat = cat(3, diag([1,1,1]), diag([2,2,2]));
            expected_matrix.label = {'mat1', 'mat2'};
            testCase.verify_spinw_matrix(testCase.swobj.matrix, ...
                expected_matrix)
        end
        
        function test_add_multiple_matrices_single_call(testCase)
            mat = cat(3, eye(3), 2*eye(3));
            testCase.swobj.addmatrix('value', mat);
            expected_matrix = testCase.default_matrix;
            expected_matrix.mat = mat;
            expected_matrix.label = {'mat1', 'mat2'};
            testCase.verify_spinw_matrix(testCase.swobj.matrix, ...
                                         expected_matrix)
        end
        
        function test_add_matrix_with_same_name_overwritten(testCase)
            testCase.swobj.addmatrix('value', 1.0);
            testCase.swobj.addmatrix('value', 2.0, 'label', 'mat1');
            expected_matrix = testCase.default_matrix;
            expected_matrix.mat = 2*eye(3);
            testCase.verify_spinw_matrix(testCase.swobj.matrix, ...
                                         expected_matrix)
        end
        
        function test_add_multiple_matrices_same_label(testCase)
            % not possible to have more than two matrix.mat with same label
            % if matrices are added with addmatrix. To test this need to
            % modify spinw attribute directly
            testCase.swobj.matrix.label = {'mat', 'mat'};
            testCase.verifyError(...
                @() testCase.swobj.addmatrix('value', 1, 'label', 'mat'), ...
                'spinw:addmatrix:LabelError');
        end
        
        function test_user_supplied_label_used(testCase)
            label = 'custom';
            testCase.swobj.addmatrix('value', 1.0, 'label', label);
            expected_matrix = testCase.default_matrix;
            expected_matrix.label = {label};
            testCase.verify_spinw_matrix(testCase.swobj.matrix, ...
                                         expected_matrix)
        end

        function test_user_supplied_color_string(testCase)
            testCase.swobj.addmatrix('value', 1.0, 'color', 'blue');
            testCase.verify_spinw_matrix(testCase.swobj.matrix, ...
                                         testCase.default_matrix)
            % test color explicitly (only size and data type checked above)
            testCase.assertEqual(testCase.swobj.matrix.color', ...
                                 [0,0,int32(255)])
        end
        
    end
     
    methods (Test, TestTags = {'Symbolic'})
        function add_matrix_with_symbolic_value(testCase)
            testCase.swobj.symbolic(true);
            testCase.swobj.addmatrix('value', 1)
            expected_matrix = testCase.default_matrix;
            expected_matrix.mat = sym('mat1')*eye(3);
            testCase.verify_spinw_matrix(expected_matrix, testCase.swobj.matrix)
        end
    end

end
