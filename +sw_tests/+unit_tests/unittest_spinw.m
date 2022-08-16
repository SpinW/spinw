classdef unittest_spinw < sw_tests.unit_tests.unittest_super
    % Runs through unit tests for @spinw/spinw.m
    properties (TestParameter)
        spinw_struct_input = { ...
            % {input struct, expected output file}
            {struct('lattice', struct('angle', [pi, pi, (2*pi)/3], ...
                                     'lat_const', [2, 2, 4])), ...
             'spinw_from_struct_lat224.mat'}, ...
            {struct('lattice', struct('angle', [pi; pi; (2*pi)/3], ...
                                     'lat_const', [2; 2; 4])), ...
             'spinw_from_struct_lat224.mat'}};
        spinw_figure_input = { ...
            % {input figure, expected output file}
            {'default_structure.fig', 'spinw_default.mat'}, ...
            {'afm_chain_spec.fig', 'spinw_afm_chain.mat'}}
        spinw_file_input = { ...
            % {input cif/fst file, expected output file}
            {'YFeO3_mcphase.cif', 'spinw_from_cif_YFeO3_mcphase.mat'}, ...
            {'LaFeO3_fullprof.cif', 'spinw_from_cif_LaFeO3_fullprof.mat'}, ...
            {'BiMn2O5.fst', 'spinw_from_fst_BiMn2O5.mat'}};
    end
    methods (Test)
        function test_spinw_no_input(testCase)
            % Tests that if spinw is called with no input, a default spinw
            % object is created
            expected_spinw = testCase.load_spinw('spinw_default.mat');
            actual_spinw = spinw;
            testCase.verify_obj(actual_spinw, expected_spinw);
        end
        function test_spinw_from_struct_input(testCase, spinw_struct_input)
            % Tests that if spinw is called with struct input, the relevant
            % fields are set
            expected_spinw = testCase.load_spinw(spinw_struct_input{2});
            actual_spinw = spinw(spinw_struct_input{1});
            testCase.verify_obj(actual_spinw, expected_spinw);
        end
        function test_spinw_from_spinw_obj(testCase)
            expected_spinw = testCase.load_spinw('spinw_from_struct_lat224.mat');
            actual_spinw = spinw(expected_spinw);
            % Ensure creating from a spinw obj creates a copy
            assert(actual_spinw ~= expected_spinw);
            testCase.verify_obj(actual_spinw, expected_spinw);
        end
        function test_spinw_from_figure(testCase, spinw_figure_input)
            expected_spinw = testCase.load_spinw(spinw_figure_input{2});
            figure = testCase.load_figure(spinw_figure_input{1});
            actual_spinw = spinw(figure);
            % Ensure creating from a figure creates a copy
            assert(actual_spinw ~= expected_spinw);
            testCase.verify_obj(actual_spinw, expected_spinw);
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
            testCase.verify_obj(actual_spinw, expected_spinw);
        end
        function test_spinw_from_file_wrong_sym(testCase)
            % Test use of a symmetry not available in symmetry.dat gives
            % an appropriate error
            fname = fullfile(testCase.get_unit_test_dir(), 'cifs', 'BiMnO3.fst');
            testCase.verifyError(@() spinw(fname), 'generator:WrongInput');
        end
        function test_spinw_from_wrong_input(testCase)
            % Test creating spinw object from invalid input gives an
            % appropriate error
            testCase.verifyError(@() spinw(4), 'spinw:spinw:WrongInput');
        end
        function test_spinw_nmagext(testCase)
             swobj = testCase.load_spinw('spinw_afm_chain.mat');
             testCase.assertEqual(swobj.nmagext, 2);
        end
        function test_spinw_natom(testCase)
             swobj = testCase.load_spinw('spinw_from_cif_YFeO3_mcphase.mat');
             testCase.assertEqual(swobj.natom, 5);
        end
    end
end
