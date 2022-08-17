classdef unittest_spinw_field < sw_tests.unit_tests.unittest_super
    properties
        swobj
    end
    properties (TestParameter)
        incorrect_input = {0, [1; 1], ones(1, 4)};
    end
    methods(TestMethodSetup)
        function create_sw_model(testCase)
            testCase.swobj = spinw();
        end
    end
    methods (Test)
        function test_incorrect_shape_field_raises_error(testCase, ...
                                                         incorrect_input)
            testCase.verifyError(...
                @() field(testCase.swobj, incorrect_input), ...
                'spinw:magfield:ArraySize');
        end
        function test_default_field(testCase)
            testCase.assertEqual(field(testCase.swobj), [0 0 0]);
        end
        function test_set_field(testCase)
            field_val = [1 2 3];
            field(testCase.swobj, field_val);
            testCase.assertEqual(field(testCase.swobj), field_val);
        end
        function test_set_field_column_vector(testCase)
            field_val = [1; 2; 3];
            field(testCase.swobj, field_val);
            testCase.assertEqual(field(testCase.swobj), field_val');
        end
        function test_set_field_multiple_times(testCase)
            field(testCase.swobj, [1 2 3]);
            new_field_val = [1.5 2.5 -3.5];
            field(testCase.swobj, new_field_val);
            testCase.assertEqual(field(testCase.swobj), new_field_val);
        end
        function test_returns_spinw_obj(testCase)
            field_val = [1 2 3];
            new_swobj = field(testCase.swobj, field_val);
            % Test field is set
            testCase.assertEqual(field(testCase.swobj), field_val)
            % Also test handle to swobj is returned
            assert (testCase.swobj == new_swobj);
        end
    end
    methods (Test, TestTags = {'Symbolic'})
        function test_set_field_sym_mode_sym_input(testCase)
            testCase.swobj.symbolic(true);
            field_val = sym([pi pi/2 pi/3]);
            field(testCase.swobj, field_val);
            testCase.assertEqual(field(testCase.swobj), field_val);
        end
        function test_set_field_sym_mode_non_sym_input(testCase)
            testCase.swobj.symbolic(true);
            field_val = [pi pi/2 pi/3];
            field(testCase.swobj, field_val);
            expected_field_val = field_val*sym('B');
            testCase.assertEqual(field(testCase.swobj), expected_field_val);
        end
    end
end