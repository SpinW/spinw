classdef unittest_spinw_fitspec < sw_tests.unit_tests.unittest_super
    % Tests for fitspec - not strictly a unit test, but make sure it runs

    properties
        datafile = '';
        swobj = [];
        fitpar = struct();
    end

    methods (TestClassSetup)
        function setup_model_and_fitpars(testCase)
            % Writes out the mode data
            testCase.datafile = fullfile(tempdir, 'triAF_modes.txt');
            fid = fopen(testCase.datafile, 'w');
            fprintf(fid, '  QH   QK   QL Elim1 Elim2    I1   EN1   sig1    I2   EN2   sig2\n');
            fprintf(fid, '   1  0.2    1   0.5   5.0  22.6  2.94  0.056     0     0      0\n');
            fprintf(fid, '   1  0.4    1   0.5   5.0  65.2  2.48  0.053  13.8  3.23  0.061\n');
            fprintf(fid, '   1  0.6    1   0.5   5.0  69.5  2.52  0.058  16.3  3.15  0.054\n');
            fprintf(fid, '   1  0.8    1   0.5   5.0  22.6  2.83  0.057     0     0      0\n');
            fclose(fid);
            testCase.swobj = sw_model('triAF', 0.7);
            testCase.fitpar = struct('datapath', testCase.datafile, ...
                                     'Evect', linspace(0, 5, 51), ...
                                     'func', @(obj,p)matparser(obj,'param',p,'mat',{'J_1'},'init',1),...
                                     'xmin', 0, 'xmax', 2, 'x0', 0.7, ...
                                     'plot', false, 'hermit', true, ...
                                     'optimizer', 'simplex', ...
                                     'maxiter', 1, 'maxfunevals', 1, 'nMax', 1, 'nrun', 1);
        end
    end

    methods (TestClassTeardown)
        function del_mode_file(testCase)
            try
                delete(testCase.datafile); % Ignores errors if file already deleted
            end
        end
    end

    methods (Test)
        function test_fitspec(testCase)
            fitout = testCase.swobj.fitspec(testCase.fitpar);
            testCase.verify_val(fitout.x, 0.7, 'abs_tol', 0.25);
            testCase.verify_val(fitout.redX2, 243, 'abs_tol', 10);
        end
        function test_fitspec_twin(testCase)
            % Checks that twins are handled correctly
            swobj = copy(testCase.swobj);
            % Adds a twin with very small volume so it doesn't affect original fit
            % If twins not handled correctly, the fit will be bad.
            swobj.addtwin('axis', [1 1 1], 'phid', 54, 'vol', 0.01);
            fitout = swobj.fitspec(testCase.fitpar);
            testCase.verify_val(fitout.x, 0.7, 'abs_tol', 0.25);
            testCase.verify_val(fitout.redX2, 243, 'abs_tol', 10);
        end
    end
end
