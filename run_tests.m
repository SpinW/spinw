function result = run_tests(out_dir)
    disp(version);
    if ~exist('spinw', 'class')
        if exist('swfiles', 'dir') && exist('external', 'dir') && exist('dat_files', 'dir')
            addpath(genpath('swfiles'));
            addpath(genpath('external'));
            addpath(genpath('dat_files'));
        else
            error(['SpinW is not installed and the swfiles, external and/or ', ...
                   'dat_files drectories couldn''t be found on the current ', ...
                   'path, so the tests cannot be run.'])
        end
    end
    if nargin == 0
        out_dir = fullfile(pwd);
    end
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    import matlab.unittest.TestSuite
    import matlab.unittest.TestRunner
    import matlab.unittest.plugins.CodeCoveragePlugin
    import matlab.unittest.plugins.codecoverage.CoberturaFormat
    import matlab.unittest.selectors.HasTag
    import matlab.unittest.plugins.XMLPlugin


    suite = TestSuite.fromPackage('sw_tests', 'IncludingSubpackages', true);
    if ~sw_hassymtoolbox()
        % only run symbolic tests when the toolbox is available
        suite = suite.selectIf(~HasTag('Symbolic'));
    end
    runner = TestRunner.withTextOutput;

    % compile mex files
    sw_mex('compile', true, 'test', false, 'swtest', false);
    
    % Add coverage output
    cov_dirs = {'swfiles', 'external'};
    for i = 1:length(cov_dirs)
        reportFormat = CoberturaFormat(fullfile(out_dir, ['coverage_', cov_dirs{i}, '.xml']));
        coverage_plugin = CodeCoveragePlugin.forFolder(cov_dirs{i}, ...
                                                       'Producing', reportFormat, ...
                                                       'IncludingSubfolders', true);
        runner.addPlugin(coverage_plugin);
        if verLessThan('matlab', '9.12') % Can add cov for multiple folders only from R2022a
            break;
        end
    end

    % Add JUnit output - unique name so they are not overwritten on CI
    junit_fname = ['junit_report_', computer, version('-release'), '.xml'];
    junit_plugin = XMLPlugin.producingJUnitFormat(junit_fname);
    runner.addPlugin(junit_plugin)

    result = runner.run(suite)
    if(any(arrayfun(@(x) x.Failed, result)))
        error('Test failed');
    end
end
