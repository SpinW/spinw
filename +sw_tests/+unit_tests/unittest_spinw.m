classdef unittest_spinw < matlab.mock.TestCase
    % Runs through unit tests for @spinw/spinw.m
    methods
        function obj = load_spinw(testCase, filename)
            obj = load(fullfile('.', 'test_data', 'unit_tests', 'spinw', filename));
            obj = obj.data;
        end
        function verify_obj(testCase, expected_obj, actual_obj, relTol, absTol)
            if nargin < 5
                absTol = 0;
            end
            if nargin < 4
                relTol = 1e-10;
            end
            all_fieldnames = union(fieldnames(expected_obj), fieldnames(actual_obj));
            for i=1:length(all_fieldnames)
                field = all_fieldnames(i);
                if strcmp(field{:}, "cache")
                    continue;
                end
                disp(field{:});
                expected_value = expected_obj.(field{:});
                actual_value = actual_obj.(field{:});
                if isstruct(expected_value)
                    testCase.verify_obj(expected_value, actual_value, relTol, absTol);
                else
                    testCase.verify_val(expected_value, actual_value, relTol, absTol);
                end
            end
        end
        function verify_val(testCase, expected_val, actual_val, relTol, absTol)
            if nargin < 5
                absTol = 0;
            end
            if nargin < 4
                relTol = 1e-10;
            end
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.RelativeTolerance
            import matlab.unittest.constraints.AbsoluteTolerance
            theseBounds = RelativeTolerance(relTol) | AbsoluteTolerance(absTol);
            testCase.verifyThat(actual_val, IsEqualTo(expected_val, 'Within', theseBounds));
        end
    end
    methods (Test)
        function test_no_input(testCase)
            % Tests that if spinw is called with no input, a default spinw
            % object is created
            expected_spinw = testCase.load_spinw('spinw_default.mat');
            actual_spinw = spinw;
            testCase.verify_obj(expected_spinw, actual_spinw);
        end
        function test_struct_input(testCase)
            % Tests that if spinw is called with struct input, the relevant
            % fields are set
            angle = [pi, pi, (2*pi)/3];
            lat_const = [2, 2, 4];
            expected_spinw = testCase.load_spinw('spinw_default.mat');
            expected_spinw.lattice.angle = angle;
            expected_spinw.lattice.lat_const = lat_const;

            actual_spinw = spinw(struct('lattice', struct('angle', angle, ...
                                                          'lat_const', lat_const)));
            testCase.verify_obj(expected_spinw, actual_spinw);
        end
    end

end
