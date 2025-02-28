function uiopen(type,direct)
% Overloaded UIOPEN to enable drag&drop for .cif files via SpinW and to
% show images in Matlab.
%
%You can drag&drop .cif files into the Matlab Command Window and
%they will be plotted automatically. For any other input, the function
%calls the original Matlab function, see it's help below.
%
%UIOPEN Present file selection dialog with appropriate file filters.
%
%   UIOPEN presents a file selection dialog.  The user can either choose a
%   file to open or click cancel.  No further action is taken if the user
%   clicks on cancel.  Otherwise the OPEN command is evaluated in the base
%   workspace with the user specified filename.
%
%   These are the file filters presented using UIOPEN.
%
%   1st input argument   Filter List
%   <no input args>      *.m, *.fig, *.mat,
%                        *.mdl, *.slx  (if Simulink is installed),
%                        *.cdr         (if Stateflow is installed),
%                        *.rtw, *.tmf, *.tlc, *.c, *.h, *.ads, *.adb
%                                      (if Simulink Coder is installed),
%                        *.*
%   MATLAB               *.m, *.fig, *.*
%   LOAD                 *.mat, *.*
%   FIGURE               *.fig, *.*
%   SIMULINK             *.mdl, *.slx, *.*
%   EDITOR               *.m, *.mdl, *.cdr, *.rtw, *.tmf, *.tlc, *.c, *.h, *.ads, *.adb, *.*
%
%   If the first input argument is unrecognized, it is treated as a file
%   filter and passed directly to the UIGETFILE command.
%
%   If the second input argument is true, the first input argument is
%   treated as a filename.
%
%   Examples:
%       uiopen % displays the dialog with the file filter set to all MATLAB
%              %files.
%
%       uiopen('matlab') %displays the dialog with the file
%                         %filter set to all MATLAB files.
%
%       uiopen('load') %displays the dialog with the file filter set to
%                      %MAT-files (*.mat).
%
%       uiopen('figure') %displays the dialog with the file filter set to
%                        %figure files (*.fig).
%
%       uiopen('simulink') %displays the dialog with the file filter set to
%                          %model files (*.mdl,*.slx).
%
%       uiopen('editor') %displays the dialog with the file filter set to
%                        %"All MATLAB files". This filters out binary
%                        %formats: .mat, .fig, .slx.
%                        %All files are opened in the MATLAB Editor.
%
%   See also UIGETFILE, UIPUTFILE, OPEN, UIIMPORT.

persistent uiroot
if isempty(uiroot)
    uiroot = [matlabroot filesep 'toolbox' filesep 'matlab' filesep 'uitools'];
    vervec = cellfun(@(x) sscanf(x, '%d'), split(version, '.'));
    if (vervec(1) * 10 + vervec(2)) >= 242
        uiroot = [uiroot filesep 'uitools'];
    end
end

if numel(type)>3 && strcmpi(type(end+(-3:0)),'.cif') && direct
    model = spinw(type);
    if ~isempty(model.matom.S)
        model.quickham([1 1 1]);
        a0 = symbol('a');
        d = [couplingtable(model,1).bondv(end,1) couplingtable(model,2).bondv(end,1) couplingtable(model,3).bondv(end,1)];
        fprintf(['The length of the first 3 inequivalent magnetic bonds are: %5.3f ' a0 ', %5.3f ' a0 ', %5.3f ' a0 '.\n'],d);
    else
        fprintf('No magnetic atom is found in the structure.\n')
    end
    plot(model);
    % add the model to base workspace
    assignin('base', 'cifmodel', model);
    fprintf('The imported SpinW object is stored in the ''cifmodel'' variable.\n');
elseif numel(type)>3 && any(strcmpi(type(end+(-3:0)),{'.png' '.jpg'})) && direct
    fName = regexp(type,[filesep '([\w\.]+?)$'],'tokens');
    figure('ToolBar','none','MenuBar','none','Name',fName{1}{1});
    warning('off','images:initSize:adjustingMag')
    imshow(type);
    warning('on','images:initSize:adjustingMag')
    %assignin('base', 'img', imread(type));
    %fprintf('The imported image is stored in the ''img'' variable.\n');
else
    pwd0 = pwd;
    cd(uiroot)
    feval('uiopen',type,direct);
    cd(pwd0)
end