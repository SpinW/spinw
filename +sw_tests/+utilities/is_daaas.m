function out = is_daaas()
    if ispc || ismac
        out = false;
    else
        [~, hostname] = system('hostname');
        % DaaaS systems have hostname of the form 'host-NNN-NNN-NNN-NNN'
        % where the number is the (internal) IP address
        out = ~isempty(regexp(hostname, 'host-[0-9\-]*', 'match'));
    end
end