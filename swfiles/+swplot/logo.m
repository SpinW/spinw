function logo(varargin)
% creates the SpinW logo
% 
% ### Syntax
% 
% `swplot.logo`
% 
% `swplot.logo(fName)`
%
% ### Description
% 
% `swplot.logo` creates and displays the SpinW logo with credentials. The
% logo is using an honest colormap [cm_inferno] and removed the tiles of
% the sine wave to symbolize the increase of quality of code (measured as a
% number of eliminated for loops) :D The logo is used for SpinW 3.0.
% Colormap is expected to change for every major version jump.
% 
% ### Examples
%
% This is the logo:
%
% ```
% >>swplot.logo
% >>snapnow
% ```
%
% ### Input Arguments
% 
% `fName`
% : File name to save the logo. Optional, if not given the logo
%   will be shown in a new figure.
% 
% ### See Also
% 
% [spinw]
%

% determine the version of SpinW
ver = sw_version;

if ~isempty(ver.Version)
    verStr = [' ' ver.Version ' ' 'R' ver.Release];
else
    verStr = [' ' 'R' ver.Release];
end

hFigure = figure('menubar','none','toolbar','none','name',...
    ['About SpinW' verStr],'NumberTitle','off','resize','off','Visible','off');
hold on

% plot wave
ds = 1.8;
%x  = linspace(0,2/ds,41);
%y  = linspace(0,0.4,11);
x  = linspace(0,2/ds,401);
y  = linspace(0,0.4,101);

[xx,yy] = ndgrid(x,y);

z = 0.2*sin(ds*xx*2*pi+pi/2);

surf(xx,yy,z,'edgecolor','none');
colormap(cm_inferno)
caxis([-0.4 0.25])
axis equal
axis off
hold on
% plot3(xx(:,1)',yy(:,1)',z(:,1)','k-')
% plot3(xx(:,end)',yy(:,end)',z(:,end)','k-')
% plot3(xx(1,:)',yy(1,:)',z(1,:)','k-')
% plot3(xx(end,:)',yy(end,:)',z(end,:)','k-')
% plot3((xx(1,:)'+xx(end,:)')/2,yy(end,:)',z(end,:)','k-')

ax = axis;
axis(ax*2)

% set window size
wSize = get(gcf,'position');
set(gcf,'position',[wSize(1:2) 560 420*0.8])
view(3)
set(gcf,'color','w')

% add name
text(-2.3,-1.3,2.35,'Spin','fontsize',110,'fontname','Krungthep')

set(hFigure,'Visible','on');

% save logo
if nargin > 0 && ~isempty(varargin{end})
    print(varargin{end},'-dpng','-r300')
    close(hFigure)
else
    % add text
    ver0 = evalc('sw_version');
    cLoc = find(ver0==',');
    if ~isempty(cLoc)
        cLoc = cLoc(1);
        ver0 = sprintf([ver0(1:cLoc) '\n' ver0((cLoc+2):end)]);
    end
    
    txt0 = sprintf([ver0 '\nWritten by:\n  S' char(225) 'ndor T' char(243) ...
        'th\n  sandor.toth@psi.ch\n  Paul Scherrer Institut\n\nContributed:\n  Simon Ward\n  Mechthild Enderle\n  Bj' char(246) 'rn F' ...
        char(229) 'k\n  Duc Manh Le\n  Rebecca Fair\n\n' ...
        'Including contributions from many authors through\nMatlab File Exchange:\n'...
        '  fminsearchnd\n  fireprint\n\n'...
        'GNU General Public License\n'...
        'You may copy, distribute and modify\nthe software as long as you '...
        'track changes/dates\nin source files. Any modifications to or '...
        'software\nincluding (via compiler) GPL-licensed code must\nalso be '...
        'made available under the GPL along with\nbuild & install instructions.']);
    tt = text(0.65,0.5,txt0,'units','normal');
    set(tt,'FontSize',9.2);
    
    WinOnTop(hFigure, true);
end

end

function WasOnTop = WinOnTop(FigureHandle, IsOnTop)
%WINONTOP allows to trigger figure's "Always On Top" state
% INPUT ARGUMENTS:
% * FigureHandle - Matlab's figure handle, scalar
% * IsOnTop      - logical scalar or empty array
%
% USAGE:
% * WinOnTop( hfigure, bool );
% * WinOnTop( hfigure );            - equal to WinOnTop( hfigure,true);
% * WinOnTop();                     - equal to WinOnTop( gcf, true);
% * WasOnTop = WinOnTop(...);       - returns boolean value "if figure WAS on top"
% * IsOnTop  = WinOnTop(hfigure,[]) - gets "if figure is on top" property
%
% LIMITATIONS:
% * java enabled
% * figure must be visible
% * figure's "WindowStyle" should be "normal"
%
%
% Written by Igor
% i3v@mail.ru
%
% 16 June 2013 - Initial version
% 27 June 2013 - removed custom "ishandle_scalar" function call
%

% Parse Inputs

if ~exist('FigureHandle','var')
    FigureHandle = gcf;
end

assert(...
    isscalar( FigureHandle ) && ishandle( FigureHandle ) &&  strcmp(get(FigureHandle,'Type'),'figure'),...
    'WinOnTop:Bad_FigureHandle_input',...
    '%s','Provided FigureHandle input is not a figure handle'...
    );

assert(...
    strcmp('on',get(FigureHandle,'Visible')),...
    'WinOnTop:FigInisible',...
    '%s','Figure Must be Visible'...
    );

assert(...
    strcmp('normal',get(FigureHandle,'WindowStyle')),...
    'WinOnTop:FigWrongWindowStyle',...
    '%s','WindowStyle Must be Normal'...
    );

if ~exist('IsOnTop','var');IsOnTop=true;end
assert(...
    islogical( IsOnTop ) &&  isscalar( IsOnTop) || isempty( IsOnTop ), ...
    'WinOnTop:Bad_IsOnTop_input',...
    '%s','Provided IsOnTop input is neither boolean, nor empty'...
    );

% Pre-checks
error(javachk('swing',mfilename)) % Swing components must be available.

% Action
% Flush the Event Queue of Graphic Objects and Update the Figure Window.
drawnow expose

warnStruct=warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
jFrame = get(handle(FigureHandle),'JavaFrame');
warning(warnStruct.state,'MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');

drawnow
v=ver('Matlab');

if str2double(v.Version) >= 8.4
    jFrame_fHGxClient = jFrame.fHG2Client;
else
    jFrame_fHGxClient = jFrame.fHG1Client;
end

WasOnTop = jFrame_fHGxClient.getWindow.isAlwaysOnTop;

if ~isempty(IsOnTop)
    jFrame_fHGxClient.getWindow.setAlwaysOnTop(IsOnTop);
end

end
