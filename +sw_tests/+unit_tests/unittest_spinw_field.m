classdef unittest_spinw_field < sw_tests.unit_tests.unittest_super
    methods (Test)
        function test_incorrect_shape_field_raises_error(testCase)
            swobj = spinw();
            testCase.verifyError(...
                @() field(swobj, [0; 0]), ...
                'spinw:magfield:ArraySize');
        end
        function test_default_field(testCase)
            swobj = spinw();
            testCase.assertEqual(field(swobj), [0 0 0]);
        end
        function test_set_field(testCase)
            swobj = spinw();
            field_val = [1 2 3];
            field(swobj, field_val);
            testCase.assertEqual(field(swobj), field_val);
        end
        function test_set_field_column_vector(testCase)
            swobj = spinw();
            field_val = [1; 2; 3];
            field(swobj, field_val);
            testCase.assertEqual(field(swobj), field_val');
        end
        function test_set_field_multiple_times(testCase)
            swobj = spinw();
            field(swobj, [1 2 3]);
            new_field_val = [1.5 2.5 -3.5];
            field(swobj, new_field_val);
            testCase.assertEqual(field(swobj), new_field_val);
        end
        function test_returns_spinw_obj(testCase)
            swobj = spinw();
            field_val = [1 2 3];
            new_swobj = field(swobj, field_val);
            % Test field is set
            testCase.assertEqual(field(swobj), field_val)
            % Also test swobj is returned
            testCase.verify_obj(swobj, new_swobj);
        end
    end
    methods (Test, TestTags = {'Symbolic'})
        function test_set_field_sym_mode_sym_input(testCase)
            swobj = spinw();
            swobj.symbolic(true);
            field_val = sym([pi pi/2 pi/3]);
            field(swobj, field_val);
            testCase.assertEqual(field(swobj), field_val);
        end
        function test_set_field_sym_mode_non_sym_input(testCase)
            swobj = spinw();
            swobj.symbolic(true);
            field_val = [pi pi/2 pi/3];
            field(swobj, field_val);
            expected_field_val = field_val*sym('B');
            testCase.assertEqual(field(swobj), expected_field_val);
        end
    end
end