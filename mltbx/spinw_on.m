if exist('spinw', 'class') ~= 8
    currdir = fileparts(mfilename('fullpath'));
    addpath(genpath(fullfile(currdir, 'swfiles')));
    addpath(genpath(fullfile(currdir, 'external')));
    addpath(genpath(fullfile(currdir, 'dat_files')));
end
