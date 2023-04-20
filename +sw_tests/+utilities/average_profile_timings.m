function average_profile_timings(save_dir)
    % function to average tic/toc timings of performance tests
    %
    % file structure produced by profile_spinwave looks like:
    % save_dir
    % |- TestName_param_1_value1_param2_value2
    %    |- 01
    %    |- 02
    %       |- tictoc_times_profile_0.txt
    %       |- tictoc_times_profile_1.txt
    arguments
        save_dir string
    end
    % get list of sub-directories (TestName_param_1_value1_param2_value2)
    test_dirs = dir(save_dir);
    test_dirs = test_dirs([test_dirs.isdir]);
    % get list unique test names (want one file per test per profile
    % setting
    test_names = {};
    for idir = 1:numel(test_dirs)
        if contains(test_dirs(idir).name, '_')
            parts = split(test_dirs(idir).name,'_');
            name = parts{1};
            if ~any(strcmp(test_names, name))
                test_names = [test_names, name];
            end
        end
    end

    for itest = 1:numel(test_names)
        % get directory names that contain test name (e.g. FMChain)
        tests = dir(fullfile(save_dir, ...
                             sprintf('*%s*', test_names{itest})));
        tests = tests([tests.isdir]);
        max_test_name = max(cellfun(@numel, {tests.name}));
        fmt_str = ['%-', num2str(max_test_name + 2, '%.0f'), 's'];
        for do_profile = 0:1
            lines = {};
            for idir = 1:numel(tests)
                % search subdirectories for tictoc files
                files = dir(fullfile(tests(idir).folder, ...
                                     tests(idir).name, ...
                                     "*", sprintf('tictoc*%.0f.txt', ...
                                                  do_profile)));
                if  ~isempty(files)
                    times = [];
                    for ifile = 1:numel(files)
                        contents = importdata(fullfile(files(ifile).folder, ...
                                                       files(ifile).name));
                        times = [times contents.data];
                    end
                    col_names = contents.textdata(2:end,1);
                    % col = avg (std) with 1 col for each func measured
                    time_str = sprintf('%.4e(%.4e)\t', ...
                                       [mean(times,2) std(times,0,2)]');
                    % add test dir name to beginning of each line
                    line = sprintf([fmt_str, '\t%s'], tests(idir).name, ...
                                   time_str);
                    lines = [lines line];
                end
            end
            % add header line at beginning (now we know how many funcs
            % measured
            lines = [['# Timings given as avg(stdev) in seconds ' ...
                      'for each function (see col. headers)'], ...
                     sprintf([fmt_str, ...
                              repmat('\t%-22s', 1, numel(col_names))], ...
                             'Test Dir.', col_names{:}), ...
                     lines];
            % write to file
            save_file = fullfile(save_dir, ...
                                 sprintf('%s_profile_%.0f.txt', ...
                                         test_names{itest}, do_profile));
            % writelines(lines, save_file);
            fid = fopen(save_file, 'w');
            for line = lines
                fprintf(fid, [line{1}, '\n']);
            end
            fclose(fid);
        end
    end
end

