function outStr = sw_version()
% returns the installed version of SpinW
% 
% ### Syntax
% 
% `sw_version`
%
% `ver = sw_version`
% 
% ### Description
% 
% `sw_version` returns the installed version of SpinW and the matlab
% version. This version number is identical to the tag of the [GitHub SpinW
% repository](https://github.com/tsdev/spinw).
%
% `ver = sw_version` returns the version information in a struct, that
% contains the program name, version, author, contact, release number,
% release date and license.
%

% Take into account deployed installs
if isdeployed
    outStr = struct;
    return
end

% read file header from sw_version.m file
fid = fopen('sw_version.m');

% first line
fgets(fid);
% skip comments
fLine = strtrim(fgets(fid));
while numel(fLine)>0 && strcmp(fLine(1),'%')
    fLine = strtrim(fgets(fid));
end
% skip empty lines
while numel(fLine)==0 || ~strcmp(fLine(1),'%')
    fLine = strtrim(fgets(fid));
end
% read version information
verLine = {fLine};
while numel(verLine{end})>0 && strcmp(verLine{end}(1),'%')
    verLine{end+1} = strtrim(fgets(fid)); %#ok<*AGROW>
end
verLine = verLine(1:end-1);
fclose(fid);

% read strings enclosed with $ signs
partStr = {};
for ii = 1:numel(verLine)
    verSel = verLine{ii};
    while sum(verSel=='$') > 1
        [~, verSel] = strtok(verSel,'$'); %#ok<*STTOK>
        [partStr{end+1}, verSel] = strtok(verSel,'$');
    end
end

nField = numel(partStr);
fieldName = cell(1,nField);
fieldVal = cell(1,nField);

verStruct = struct;

% extract values and save them into a structure
for ii = 1:nField
    [fieldName{ii}, fieldVal{ii}] = strtok(partStr{ii},':');
    fieldVal{ii} = strtrim(fieldVal{ii}(2:end));
    verStruct.(fieldName{ii}) = fieldVal{ii};
end


if nField == 0
    aDir = pwd;
    cd(sw_rootdir);
    [~, revNum] = system('git rev-list --count HEAD');
    revNum = strtrim(revNum);
    %[~, revNum] = system('svn info |grep Revision: |cut -c11-');
    cd(aDir);
    revNum = str2double(revNum)+1e3;
end

% Matlab version & Symbolic Toolbox
if ~license('checkout','Symbolic_Toolbox')
    strSym = 'no Symbolic Math Toolbox installed';
else
    strSym = 'Symbolic Math Toolbox installed';
end


if nargout == 0
    if nField == 0
        
        if any(revNum)
            fprintf('This version of SpinW (rev. num. %d) is not released yet!\n',revNum);
        else
            fprintf('This version of SpinW is not released yet!\n');
        end
    else
        disp([verStruct.Name verStruct.Version ' (rev ' num2str(verStruct.Release) ')']);
        onlineRev = sw_update;
        if onlineRev > str2num(verStruct.Release) %#ok<ST2NM>
            disp(['Newer version of SpinW is available online (rev. num. ' num2str(onlineRev) '), use the sw_update() function to download it!']);
        else
            disp('You have the latest version of SpinW!')
        end
    end
    fprintf(['MATLAB version: ' version ', ' strSym '\n']);
    
else
    ver0 = struct;
    ver0.Name     = 'SpinW';
    ver0.Version  = '';
    ver0.Release  = '';
    ver0.Date     = datestr(now,'dd-mmm-yyyy');
    ver0.Author   = 'S. Tóth and S. Ward';
    ver0.Contact  = 'spinw4@gmail.com, @spinw4 on Twitter';
    ver0.License  = 'GNU GENERAL PUBLIC LICENSE';

    if nField == 0
        if any(revNum)
            ver0.Release = num2str(revNum);
        end
        outStr = ver0;
    else
        if isempty(fieldnames(verStruct))
            outStr = ver0;
        else
            outStr = verStruct;
        end
    end
end

end