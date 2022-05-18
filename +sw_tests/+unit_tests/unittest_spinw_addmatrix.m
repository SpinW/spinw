classdef unittest_spinw_addmatrix < matlab.mock.TestCase
    % Runs through unit test for @spinw/spinwave.m

    properties
        swobj = [];
    end

    methods (TestMethodSetup)
        function setup_spinw_model(testCase)
            testCase.swobj = spinw; % default init
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
        
        function test_non_numeric_input_warns_and_not_add_matrix(testCase)
            testCase.verifyWarning(...
                @() testCase.swobj.addmatrix('value', 'str'), ...
                'spinw:addmatrix:WrongInput');
            testCase.assertTrue(isempty(testCase.swobj.matrix.mat));
        end
        
        function test_single_value_adds_diag_matrix_no_label_color(testCase)
            J = 1.0;
            testCase.swobj.addmatrix('value', J);
            % check diagonal elements are equal to J
            testCase.assertTrue(isdiag(testCase.swobj.matrix.mat));
            for elem = diag(testCase.swobj.matrix.mat)'
               testCase.assertEqual(elem, J) 
            end
            % check default label and color 3 vector (randomly assigned)
            testCase.assertEqual(testCase.swobj.matrix.label{1}, 'mat1')
            testCase.assertEqual(size(testCase.swobj.matrix.color), [3,1])
        end
        
        function test_matrix_value_added_no_modification(testCase)
            matrix = reshape(1:9, 3, 3);
            testCase.swobj.addmatrix('value', matrix);
            % check diagonal elements are equal to J
            testCase.assertEqual(testCase.swobj.matrix.mat, matrix);
        end
        
        function test_vector_value_adds_DM_matrix(testCase)
            [m1, m2, m3] = deal(1,2,3);
            testCase.swobj.addmatrix('value', [m1, m2, m3]);
            % check diagonal elements are equal to J
            testCase.assertEqual(testCase.swobj.matrix.mat, ...
                [0 m3 -m2; -m3 0 m1; m2 -m1 0]);
        end
        
        function test_add_multiple_matrices(testCase)
            nmat = 2;
            for imat = 1:nmat
                testCase.swobj.addmatrix('value', imat);
            end
            % check matrix properties have correct dimensions
            testCase.assertEqual(size(testCase.swobj.matrix.mat), [3,3,nmat])
            testCase.assertEqual(size(testCase.swobj.matrix.color), [3,nmat])
            testCase.assertEqual(size(testCase.swobj.matrix.label), [1, nmat])
            % check value and label correct
            for imat = 1:nmat
                for elem = diag(testCase.swobj.matrix.mat(:,:,imat))'
                    testCase.assertEqual(elem, imat) 
                end
                testCase.assertEqual(testCase.swobj.matrix.label{imat}, ...
                    ['mat' num2str(imat)])
            end
        end
        
        function test_user_supplied_label_used(testCase)
            testCase.swobj.addmatrix('value', 1.0, 'label', 'custom');
            testCase.assertEqual(testCase.swobj.matrix.label{1}, 'custom')
        end

        function test_user_supplied_color_string(testCase)
            testCase.swobj.addmatrix('value', 1.0, 'color', 'blue');
            testCase.assertEqual(testCase.swobj.matrix.color', [0,0,int32(255)])
        end
        
     end

end
