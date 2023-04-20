function run_performance_tests(nruns)
    if nargin == 0
        nruns = 1;
    end
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
    
    % run tests
    fpath_parts = {'+sw_tests', '+performance_tests'};
    files  = dir(fullfile(fpath_parts{:}, '*.m'));
    for fname = {files.name}
        for nrun = 1:nruns
            run(erase(strjoin([fpath_parts, fname{1}], '.'), '+'));
        end
    end

end
