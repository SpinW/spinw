function mem = sw_freemem
% calculates the available memory
% 
% ### Syntax
% 
% `mem = sw_freemem`
% 
% ### Description
% 
% `mem = sw_freemem` determines the available free memory (RAM). If the
% function cannot determine the size of the free memory, it returns zero.
% The function is compatible with Linux, macOS and Windows.
% 
% ### Output Arguments
%
% `mem`
% : Size of free memory in bytes.
%
persistent t0 m0;
if isempty(t0) || isempty(m0)
    m0 = get_free();
    t0 = tic;
else
    if isempty(t0) || toc(t0) > 10    % poll only every 10 seconds
        m0 = get_free();
        t0 = tic;
    end
end
mem = m0;

end

function mem = get_free()
mem = 0;

try %#ok<TRYNC>
    if ispc
        memStr = memory;
        mem = memStr.MemAvailableAllArrays;
    elseif ismac
        [~,memStr] = unix('vm_stat | grep free');
        mem = sscanf(memStr(14:end),'%f')*4096;
    elseif isunix
        % get the buffer/cache size
        [~, memStr] = unix('free -b | grep ''-''');
        
        if isempty(memStr)
            % there is no buffer/cache, just get the 'Mem' calues
            [~, memStr] = unix('free -b | grep ''Mem''');
            [~, mem_free] = strtok(memStr);
            mem = sscanf(mem_free,'%f');
            mem = mem(3);
        else
            [~, mem_free] = strtok(memStr(20:end));
            mem = str2double(mem_free);
        end
    end
end

end