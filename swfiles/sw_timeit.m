function sw_timeit(percent,varargin)
% timer and remaining time estimator
% 
% ### Syntax
% 
% `sw_timeit(percent, {mode},{tid},{title})`
% 
% ### Description
% 
% `sw_timeit(percent, {mode},{tid},{title})` can display remaining time of
% a calculation that is run for a fixed number of iterations. It can output
% the status both in the Command Window and in a pup up window using
% [waitbar].
% 
% ### Input Arguments
% 
% `percent`
% : Percentage of the calculation that is already done.
% 
% `mode`
% : Controls the time estimation, optional parameter:
%   * `1` Starts the time estimation.
%   * `0` Displays of the remaining time, default value.
%   * `2` The calculation finished, show a summary.
% 
% `tid`
% : Determines if the elapsed and required time for the calculation is
%   displayed. The default value is determined by the `tid` preference
%   stored in [swpref]. The following values are allowed:
%   * `0` No timing is displayed.
%   * `1` Display the timing in the Command Window.
%   * `2` Show the timing in a separat pup-up window.
%
% `title`
% : The text to show in the pup up window.
% 
% ### See Also
% 
% [waitbar]
%

global sw_time
pref = swpref;

if nargin == 0
    swhelp sw_timeit
    return
end

if nargin > 2 && ~isempty(varargin{2}) && ~ischar(varargin{2})
    fid = varargin{2};
else
    fid = pref.tid;
end

if fid == 0
    % do nothing
    return
end

if nargin>3
    title0 = varargin{3};
else
    title0 = 'sw_timeit';
end

if ~ismember(fid,[1 2])
    return
end

if nargin > 1
    start = varargin{1};
else
    start = 0;
end

switch start
    case 1
        % start the time estimation
        sw_time = tic;
        switch fid
            case 1
                fprintf([repmat(' ',[1 40]) '\n']);
            case 2
                hBar = waitbar(0,'Initializing...');
                hBar.HandleVisibility='on';
                hBar.Tag = 'sw_timeit';
                hBar.Name = title0;
                drawnow;
        end
    case 0
        % refresh the displayed time
        etime = double(toc(sw_time));
        if percent == 0
            percent = 1e-5;
        end
        rtime = (100-percent)/percent*etime;
        hou = floor(rtime/60^2);
        rtime = rtime-hou*60^2;
        min = floor(rtime/60);
        sec = floor(rtime - min*60);
        switch fid
            case 1
                fprintf([repmat('\b',[1 41]) '%6.2f%%, remained: %03d:%02d:%02d (HH:MM:SS).\n'],...
                    percent,hou,min,sec);
            case 2
                hBar = findobj('Tag','sw_timeit');
                if ~isempty(hBar)
                    waitbar(percent/100,hBar(1),sprintf('%6.2f%%, remained: %03d:%02d:%02d (HH:MM:SS)',percent,hou,min,sec))
                    drawnow;
                end
        end
        
    case  2
        % finish extimation
        etime = double(toc(sw_time));
        hou = floor(etime/60^2);
        etime = etime-hou*60^2;
        min = floor(etime/60);
        etime = etime-min*60;
        sec = floor(etime);
        %tho = floor((etime-sec)*1000);
        %fprintf('Finished in %02d:%02d:%02d.%03d (HH:MM:SS.FFF).\n',hou,min,sec,tho);
        switch fid
            case 1
                fprintf(repmat('\b',1,40+1));
                fprintf('Calculation is finished in %02d:%02d:%02d (hh:mm:ss).\n',hou,min,sec);
            case 2
                hBar = findobj('Tag','sw_timeit');
                delete(hBar);
                fprintf('Calculation is finished in %02d:%02d:%02d (hh:mm:ss).\n',hou,min,sec);
        end
end

end

function extended_waitbar()
% extended waitbar to stop execution and start debug mode

h = waitbar(0,'Progress','Name','Waitbar',...
    'CreateCancelBtn',...
    'setappdata(gcbf,''canceling'',1)');

hChild = get(h,'Children');
hCancelBtn = hChild(ismember(get(hChild,'Tag'),'TMWWaitbarCancelButton'));
pauseBtnPos = get(hCancelBtn,'Position');
pauseBtnPos(1) = pauseBtnPos(1) - pauseBtnPos(3) - 1;
pauseBtnPos(2) = pauseBtnPos(2)+1;

hPauseBtn = uicontrol(h, 'Style', 'pushbutton', 'String', 'Debug',...
    'Position', pauseBtnPos,...
    'Callback', 'dbstop in sw_timeit.m at 23');

end