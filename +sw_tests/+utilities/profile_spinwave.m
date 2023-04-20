function profile_spinwave(test_name, sw_obj, spinwave_args, egrid_args, ...
                          inst_args, do_profiles)

    if nargin < 6
        do_profiles = 1;
    end

    % generate file directory name for results
    try
        % get current commit is possible
        commit = evalc("!git rev-parse --short HEAD");
        commit = commit(1:end-1); % remove newline
    catch
        ver = sw_version();
        commit = [ver.Name ver.Release];
    end
    host_info = [computer(), '_', version('-release')];
    default_save_dir = fullfile(pwd, "profile_results", commit, ...
        host_info, test_name);
    irun = 0;
    is_dir = true;
    while is_dir
        save_dir = fullfile(default_save_dir, num2str(irun, '%02.0f'));
        is_dir = isfolder(save_dir);
        irun = irun + 1;
    end
    mkdir(save_dir);

    for do_profile = do_profiles
        % open file for tic/toc timings
        fid = fopen(fullfile(save_dir, ...
                    sprintf('tictoc_times_profile_%i.txt', do_profile)), 'w');
        c = onCleanup(@()fclose(fid));  % in case of exception file will close
        fprintf(fid, "Function\tDuration (s)\n");
    
        if do_profile
            % start profiling
            profile('clear');
            profile('on', '-memory');
        else
             profile('off'); % can be left on if user aborts prematurely
        end
        % use supercell k=0 structure
        start_time = tic;
        spec = sw_obj.spinwave(spinwave_args{:});
        fprintf(fid, 'spinwave\t%.4e\n', toc(start_time));
        if ~isempty(egrid_args)
            tic;
            spec = sw_egrid(spec, egrid_args{:});
            fprintf(fid, 'sw_egrid\t%.4e\n', toc);
            if ~isempty(inst_args)
                tic;
                sw_instrument(spec, inst_args{:});
                fprintf(fid, 'sw_instrument\t%.4e\n', toc);
            end
        end
    
        if do_profile
            % save profile results
            p = profile('info');
            profsave(p, save_dir);  % will mkdir if not exist
            % save ascii summary
            sw_tests.utilities.save_profile_results_to_txt(p, save_dir);
            profile('off');
        end
    end
end