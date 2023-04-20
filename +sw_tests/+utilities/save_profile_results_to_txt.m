function save_profile_results_to_txt(results, save_dir, fname_suffix)
    arguments
        results struct
        save_dir string
        fname_suffix string = string();  % default is empty
    end

    extract = {'FunctionName' 'NumCalls' 'TotalTime' 'TotalMemAllocated' 'TotalMemFreed' 'PeakMem'};

    ft = results.FunctionTable;
    if ~isempty(ft)
        maxTime = max([ft.TotalTime]);
    
        fn = fieldnames(ft);
        sd = setdiff(fn, extract);
        m = rmfield(ft, sd);
    
        percent = arrayfun(@(x) 100*x.TotalTime/maxTime, m, 'UniformOutput', false);
        [m.PercentageTime] = percent{:};
        sp_time = arrayfun(@(x) sum([x.Children.TotalTime]), ft);
        self_time = arrayfun(@(x,y) x-y, [ft.TotalTime]', sp_time, 'UniformOutput', false);
        [m.SelfTime] = self_time{:};
        percent = arrayfun(@(x) 100*x.SelfTime/maxTime, m, 'UniformOutput', false);
        [m.SelfPercentageTime] = percent{:};
    
        dataStr = evalc('struct2table(m)');
    
        % Remove HTML, braces and header
        dataStr = regexprep(dataStr, '<.*?>', '');
        dataStr = regexprep(dataStr, '[{}]', ' ');
        dataStr = dataStr(24:end);
    
        % make filename using function and mfilename of test
        filepath = fullfile(save_dir, 'summary.txt');
        % write profile results to file
        fh = fopen(filepath, 'w');
        fwrite(fh, dataStr);
        fclose(fh);
    end
end