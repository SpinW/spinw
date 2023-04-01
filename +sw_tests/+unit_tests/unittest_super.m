classdef unittest_super < matlab.mock.TestCase
    properties
        cleanup_warnings = {};
    end
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
        function verify_obj(testCase, actual_obj, expected_obj, varargin)
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
                    testCase.verify_obj(actual_value, expected_value, varargin{:});
                else
                    testCase.verify_val(actual_value, expected_value, ...
                        'field', field{:}, varargin{:});
                end
            end
        end
        function verify_val(testCase, actual_val, expected_val, varargin)
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
                        rel_tol = varargin{iarg + 1};
                    case 'field' 
                        field = varargin{iarg + 1};
                end
            end
            bounds = RelativeTolerance(rel_tol) | AbsoluteTolerance(abs_tol);
            testCase.verifyThat(actual_val, ...
                IsEqualTo(expected_val, 'Within', bounds), field);
        end
        function verify_spinw_matrix(testCase, actual_matrix, expected_matrix, varargin)
            % compare excl. color (which is randomly generated)
            testCase.verify_val(rmfield(actual_matrix, 'color'), ...
                rmfield(expected_matrix, 'color'), varargin{:})
            % check size and data type of color
            testCase.assertEqual(size(actual_matrix.color), ...
             [3, size(expected_matrix.mat, 3)]);
            testCase.assertTrue(isa(actual_matrix.color, ...
                class(expected_matrix.color)));
        end
        function verify_spinwave(testCase, actual_spinwave, ...
                                 expected_spinwave, varargin)
            % List of fields to test separately, only remove fields that
            % exist
            rmfields = intersect(fields(expected_spinwave), ...
                                 {'datestart', 'dateend', 'obj', 'V'});
            testCase.verify_obj(rmfield(actual_spinwave, rmfields), ...
                                rmfield(expected_spinwave, rmfields), varargin{:})
            for field = ["datestart", "dateend"]
                if isfield(expected_spinwave, field)
                    testCase.assertTrue(isa(actual_spinwave.(field), 'char'));
                end
            end
            % obj is not always in spinwave output (e.g. if fitmode ==
            % true)
            if isfield(expected_spinwave, 'obj')
                testCase.verify_obj(actual_spinwave.obj, expected_spinwave.obj);
            end
            % verify abs of V (matrix of eigenvecs) - sign doesn't matter
            % get sign by comaparing max abs value
            if isfield(expected_spinwave, 'V')
                ifinite = find(isfinite(expected_spinwave.V));
                [~, imax] = max(abs(expected_spinwave.V(ifinite)));
                imax = ifinite(imax); % get index in array incl. non-finite
                scale_sign = sign(expected_spinwave.V(imax)./actual_spinwave.V(imax));
                if ~isfinite(scale_sign)
                    % in case actual V(imax) is not finite - verify should fail!
                    scale_sign = 1;
                end
                testCase.verify_val(actual_spinwave.V, ...
                                    scale_sign.*expected_spinwave.V,...
                                    varargin{:});
            end
        end
        function disable_warnings(testCase, varargin)
            testCase.cleanup_warnings = [testCase.cleanup_warnings, ...
                {onCleanup(@(c) cellfun(@(c) warning('on', c), varargin))}];
            cellfun(@(c) warning('off', c), varargin);
        end
    end
end
