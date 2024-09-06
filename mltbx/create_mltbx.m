currdir = fileparts(mfilename('fullpath'));
mltbx_dir = fullfile(currdir, 'mltbx');
if exist(mltbx_dir)
    rmdir(mltbx_dir, 's');
end
mkdir(fullfile(currdir, 'mltbx'));
copyfile(fullfile(currdir, '..', 'CITATION.cff'), fullfile(currdir, 'mltbx'));
copyfile(fullfile(currdir, '..', 'license.txt'), fullfile(currdir, 'mltbx'));
copyfile(fullfile(currdir, '..', 'swfiles'), fullfile(currdir, 'mltbx', 'swfiles'));
copyfile(fullfile(currdir, '..', 'external'), fullfile(currdir, 'mltbx', 'external'));
copyfile(fullfile(currdir, '..', 'dat_files'), fullfile(currdir, 'mltbx', 'dat_files'));
copyfile(fullfile(currdir, 'spinw_on.m'), fullfile(currdir, 'mltbx'));
matlab.addons.toolbox.packageToolbox("spinw.prj", "spinw.mltbx");
