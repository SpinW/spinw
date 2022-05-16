classdef unittest_spinw < matlab.mock.TestCase
    % Runs through unit tests for @spinw/spinw.m
    properties (TestParameter)
        spinw_struct_input = { ...
            {struct('lattice', struct('angle', [pi, pi, (2*pi)/3], ...
                                     'lat_const', [2, 2, 4])),
             'spinw_from_struct_lat224.mat'}, ...
            {struct('lattice', struct('angle', [pi; pi; (2*pi)/3], ...
                                     'lat_const', [2; 2; 4])), ...
             'spinw_from_struct_lat224.mat'}};
        spinw_figure_input = { ...
            {'default_structure.fig', 'spinw_default.mat'}, ...
            {'afm_chain_spec.fig', 'spinw_afm_chain.mat'}}
        spinw_file_input = { ...
            {'YFeO3_mcphase.cif', 'spinw_from_cif_YFeO3_mcphase.mat'}, ...
            {'LaFeO3_fullprof.cif', 'spinw_from_cif_LaFeO3_fullprof.mat'}, ...
            {'BiMn2O5.fst', 'spinw_from_fst_BiMn2O5.mat'}};
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
        function verify_obj(testCase, expected_obj, actual_obj, rel_tol, abs_tol)
            if nargin < 5
                abs_tol = 0;
            end
            if nargin < 4
                rel_tol = 1e-10;
            end
            all_fieldnames = union(fieldnames(expected_obj), fieldnames(actual_obj));
            for i=1:length(all_fieldnames)
                field = all_fieldnames(i);
                if strcmp(field{:}, "cache")
                    continue;
                end
                expected_value = expected_obj.(field{:});
                actual_value = actual_obj.(field{:});
                if isstruct(expected_value)
                    testCase.verify_obj(expected_value, actual_value, rel_tol, abs_tol);
                else
                    testCase.verify_val(expected_value, actual_value, rel_tol, abs_tol, field{:});
                end
            end
        end
        function verify_val(testCase, expected_val, actual_val, rel_tol, abs_tol, field)
            if nargin < 5
                abs_tol = 0;
            end
            if nargin < 4
                rel_tol = 1e-10;
            end
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.constraints.RelativeTolerance
            import matlab.unittest.constraints.AbsoluteTolerance
            theseBounds = RelativeTolerance(rel_tol) | AbsoluteTolerance(abs_tol);
            testCase.verifyThat(actual_val, IsEqualTo(expected_val, 'Within', theseBounds), field);
        end
    end
    methods (Test)
        function test_spinw_no_input(testCase)
            % Tests that if spinw is called with no input, a default spinw
            % object is created
            expected_spinw = testCase.load_spinw('spinw_default.mat');
            actual_spinw = spinw;
            testCase.verify_obj(expected_spinw, actual_spinw);
        end
        function test_spinw_from_struct_input(testCase, spinw_struct_input)
            % Tests that if spinw is called with struct input, the relevant
            % fields are set
            expected_spinw = testCase.load_spinw(spinw_struct_input{2});
            actual_spinw = spinw(spinw_struct_input{1});
            testCase.verify_obj(expected_spinw, actual_spinw);
        end
        function test_spinw_from_spinw_obj(testCase)
            expected_spinw = testCase.load_spinw('spinw_from_struct_lat224.mat');
            actual_spinw = spinw(expected_spinw);
            % Ensure creating from a spinw obj creates a copy
            assert(actual_spinw ~= expected_spinw);
            testCase.verify_obj(expected_spinw, actual_spinw);
        end
        function test_spinw_from_figure(testCase, spinw_figure_input)
            expected_spinw = testCase.load_spinw(spinw_figure_input{2});
            figure = testCase.load_figure(spinw_figure_input{1});
            actual_spinw = spinw(figure);
            % Ensure creating from a figure creates a copy
            assert(actual_spinw ~= expected_spinw);
            testCase.verify_obj(expected_spinw, actual_spinw);
            close(figure);
        end
        function test_spinw_from_incorrect_figure(testCase)
            fig = figure('visible', 'off');
            testCase.verifyError(@() spinw(fig), 'spinw:spinw:WrongInput');
            close(fig);
        end
        function test_spinw_from_file(testCase, spinw_file_input)
            fname = fullfile(testCase.get_unit_test_dir(), 'cifs', spinw_file_input{1});
            expected_spinw = testCase.load_spinw(spinw_file_input{2});
            actual_spinw = spinw(fname);
            testCase.verify_obj(expected_spinw, actual_spinw);
        end
        function test_spinw_from_file_wrong_sym(testCase)
            % Test use of a symmetry not available in symmetry.dat gives
            % an appropriate error
            fname = fullfile(testCase.get_unit_test_dir(), 'cifs', 'BiMnO3.fst');
            testCase.verifyError(@() spinw(fname), 'generator:WrongInput');
        end
    end

end
