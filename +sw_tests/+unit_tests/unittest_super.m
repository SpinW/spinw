classdef unittest_super < matlab.mock.TestCase
    methods (Static)
        function udir = get_unit_test_dir()
            udir = fullfile('.', 'test_data', 'unit_tests');
        end
    end
    methods
        function obj = load_spinw(testCase, filename)
            obj = load(fullfile(testCase.get_unit_test_dir(), 'spinw', filename));
            obj = obj.data;
        end
        function obj = load_figure(testCase, filename)
            path = fullfile(testCase.get_unit_test_dir(), 'Figure', filename);
            obj = openfig(path, 'invisible');
        end
        function verify_obj(testCase, expected_obj, actual_obj, varargin)
            testCase.assertClass(actual_obj, class(expected_obj));
            all_fieldnames = fieldnames(expected_obj);
            if isa(expected_obj, 'struct')
                all_fieldnames = union(all_fieldnames, fieldnames(actual_obj));
            end
            for i=1:length(all_fieldnames)
                field = all_fieldnames(i);
                if strcmp(field{:}, "cache")
                    continue;
                end
                expected_value = expected_obj.(field{:});
                actual_value = actual_obj.(field{:});
                if isstruct(expected_value)
                    testCase.verify_obj(expected_value, actual_value, varargin{:});
                else
                    testCase.verify_val(expected_value, actual_value, ...
                        'field', field{:}, varargin{:});
                end
            end
        end
        function verify_val(testCase, expected_val, actual_val, varargin)
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.RelativeTolerance
            import matlab.unittest.constraints.AbsoluteTolerance
            
            field = "";
            abs_tol = 10*eps; % abs_tol comparsion is default
            rel_tol = 0;
            for iarg = 1:2:numel(varargin)
                switch varargin{iarg}
                    case 'abs_tol'
                        abs_tol = varargin{iarg + 1};
                    case 'rel_tol'
                        if expected_val > 0
                            rel_tol = varargin{iarg + 1};
                        end
                    case 'field' 
                        field = varargin{iarg + 1};
                end
            end
            bounds = RelativeTolerance(rel_tol) | AbsoluteTolerance(abs_tol);
            testCase.verifyThat(actual_val, ...
                IsEqualTo(expected_val, 'Within', bounds), field);
        end
    end
end
