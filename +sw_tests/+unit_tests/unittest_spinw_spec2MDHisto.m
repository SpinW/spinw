classdef unittest_spinw_spec2MDHisto < sw_tests.unit_tests.unittest_super
    properties
        swModel = [];
        tmpdir = '';
        testfilename = '';
        nsteps = {100};
    end
    
    properties (TestParameter)
        testpars = struct(...
            'test_1_0_0', struct('q0', [-3 0 0], 'qmax', [1 0 0], 'proj', [[1 0 0]' [0 1 0]' [0 0 1]'], 'nxs', 'test100mdh.nxs'), ...
            'test_1_1_0', struct('q0', [-3 -3 0], 'qmax', [1 1 0], 'proj', [[1 1 0]' [1 -1 0]' [0 0 1]'], 'nxs', 'test110mdh.nxs'), ...
            'test_1_1_1', struct('q0', [0 0 0], 'qmax', [1 1 1], 'proj', [[1 1 1]' [1 -1 0]' [1 1 -2]'], 'nxs', 'test111mdh.nxs'), ...
            'test_1_1_2', struct('q0', [0 0 2], 'qmax', [1 1 2], 'proj', [[1 1 0]' [1 -1 0]' [0 0 1]'], 'nxs', 'test112mdh.nxs'), ...
            'test_1_1_2_2', struct('q0', [0 0 2], 'qmax', [1 1 2], 'proj', [[1 -1 0]' [1 1 0]' [0 0 1]'], 'nxs', 'test112_2mdh.nxs'), ...
            'test_2_2_2', struct('q0', [2 2 2], 'qmax', [3 3 2], 'proj', [[1 1 0]' [1 -1 0]' [0 0 1]'], 'nxs', 'test222mdh.nxs'));
    end

    methods (TestClassSetup)
        function setup_model(testCase)
            testCase.swModel = sw_model('triAF', 1);
        end
        function setup_tempdir(testCase)
            testCase.tmpdir = tempdir;
        end
    end

    methods (TestMethodTeardown)
        function remove_tmpdir(testCase)
            delete(testCase.testfilename);
        end
    end

    methods (Test)
        function test_qdirs(testCase, testpars)
            q0 = testpars.q0;
            qmax = testpars.qmax;
            proj = testpars.proj;
            testCase.disable_warnings('spinw:spinwave:NonPosDefHamiltonian');
            spec = sw_egrid(spinwave(testCase.swModel, {q0 qmax testCase.nsteps{1}}));
            % dproj = [(qmax-q0)/testCase.nsteps{1}, 1e-6, 1e-6];
            dproj = [1, 1e-6, 1e-6];
            testCase.testfilename = fullfile(testCase.tmpdir, testpars.nxs);
            sw_spec2MDHisto(spec, proj, dproj, testCase.testfilename);
        end

        function test_non_ortho(testCase)
            q0 = [0 0 2];
            qdir = [1 1 0];
            spec = sw_egrid(spinwave(testCase.swModel, {q0 q0+qdir testCase.nsteps{1}}));
            proj = [qdir(:) [1 0 0]' [0 0 1]'];
            dproj = [1e-6, norm((qdir-q0))/testCase.nsteps{1}, 1e-6];
            verifyError(testCase,@() sw_spec2MDHisto(spec, proj, dproj, 'tmp/test_blank.nxs'), "read_struct:nonorthogonal")
        end
    end
end