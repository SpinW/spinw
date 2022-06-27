classdef unittest_spinw_addtwin < sw_tests.unit_tests.unittest_super

    properties
        swobj = [];
    end

    methods (TestMethodSetup)
        function setup_spinw_model(testCase)
            testCase.swobj = spinw();
        end
    end

    methods (Test)
        function test_no_input_calls_help(testCase)
            % Tests that if call addmatrix with no input, it calls the help
            help_function = sw_tests.utilities.mock_function('swhelp');
            testCase.swobj.addtwin();
            testCase.assertEqual(help_function.n_calls, 1);
            testCase.assertEqual(help_function.arguments, ...
                {{'spinw.addtwin'}});
        end
        
        function test_addtwin_with_rotc_not_valid_rotation(testCase)
            testCase.verifyError(...
                @() testCase.swobj.addtwin('rotC', ones(3)), ...
                'spinw:addtwin:WrongInput');
        end
        
        function test_addtwin_with_rotc_invalid_size(testCase)
            testCase.verifyError(...
                @() testCase.swobj.addtwin('rotC', ones(4) ), ...
                'sw_readparam:ParameterSizeMismatch');
        end
        
        function test_addtwin_with_only_axis_defaults_identity_rotc(testCase)
            testCase.swobj.addtwin('axis', [0 0 1])
            testCase.assertEqual(testCase.swobj.twin.vol, [1, 1])
            testCase.assertEqual(testCase.swobj.twin.rotc, ...
                cat(3, eye(3), eye(3)))
        end
        
        function test_addtwin_with_axis_custom_vol_ratio(testCase)
            testCase.swobj.addtwin('axis', [0 0 1], 'vol', 0.5)
            testCase.assertEqual(testCase.swobj.twin.vol, [1, 0.5])
        end
        
        function test_addtwin_overwrite(testCase)
            testCase.swobj.addtwin('axis', [0,0,1], 'phid', 90, ...
                'vol', 0.5, 'overwrite', true)
            testCase.assertEqual(testCase.swobj.twin.vol, 0.5)
            testCase.verify_val(testCase.swobj.twin.rotc, ...
                [0 -1 0; 1 0 0; 0 0 1])
        end
        
        function test_addtwin_with_phid(testCase)
            import matlab.unittest.constraints.IsEqualTo
            testCase.swobj.addtwin('axis', [0 0 1], 'phid', [90])
            testCase.assertEqual(testCase.swobj.twin.vol, [1, 1])          
            testCase.verify_val(testCase.swobj.twin.rotc(:,:,2), ...
                [0 -1 0; 1 0 0; 0 0 1])
        end
          
        function test_addtwin_with_multiple_phid(testCase)
            testCase.swobj.addtwin('axis', [0 0 1], 'phid', [90, 180])
            testCase.assertEqual(testCase.swobj.twin.vol, [1, 1, 1])          
            testCase.verify_val(testCase.swobj.twin.rotc(:,:,2), ...
                [0 -1 0; 1 0 0; 0 0 1])
            testCase.verify_val(testCase.swobj.twin.rotc(:,:,3), ...
               [-1 0 0; 0 -1 0; 0 0 1])
            
        end
        
        function test_addtwin_with_multiple_rotc(testCase)
            rotc = cat(3, [0 -1 0; 1 0 0; 0 0 1], [-1 0 0; 0 -1 0; 0 0 1]);
            testCase.swobj.addtwin('rotc', rotc)
            testCase.assertEqual(testCase.swobj.twin.vol, [1, 1, 1])
            testCase.assertEqual(testCase.swobj.twin.rotc(:,:,2:end), rotc)  
        end
        
        function test_addtwin_overwrite_vol_ratio(testCase)
            testCase.swobj.addtwin('axis', [0,0,1], 'vol', 0.5, ...
                'overwrite', true)
            testCase.assertEqual(testCase.swobj.twin.vol, 0.5)
        end
        
    end
end
