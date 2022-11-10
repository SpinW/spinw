function [fHandle0, pHandle0] = sw_plotspec(spectra, varargin)
% plots spectrum
% 
% ### Syntax
% 
% `[fhandle, phandle] = sw_plotspec(spectra,Name,Value)`
% 
% ### Description
% 
% `[fhandle, phandle] = sw_plotspec(spectra,Name,Value)` plots excitation
% spectrum that is calculated either by [spinw.spinwave] or
% [spinw.powspec]. It can plot dispersion or intensities as line plots or
% the energy binned spectrum as a color plot. The color plots uses
% [cm_inferno] as a default colormap. To change the default colormap use
% the `swpref.setpref('colormap',@my_colomap)` command. The function is
% able to plot the spectrum if it is calculated along a path in the
% Brillouin-zone and display the labels of the high symmetry Brillouon-zone
% points.
% 
% ### Name-Value Pair Arguments
% 
% `'mode'`
% : Choose the type of plot using the following strings:
%   * `'disp'`  Plot dispersion as line plot.
%   * `'int'`   PLot intensity of each mode as line plot.
%   * `'color'` Color plot of energy binned spectrum.
%   * `'auto'`  Auto plot mode that tries to determine the best
%               parameteres, default.
% 
% `'imag'`
% : If `true` also the imaginary part of the dispersion
%   and the correlation function values will be shown as red lines on top
%   the real values. For color plot if `true` only the imaginary part of
%   the binned data will be shown. Default value is `false`.
% 
% `'aHandle'`
% : Handle of the axis object which will show the plot. If undefined the
%   active axis will be used, see [gca].
% 
% `'colorbar'`
% : Plot colorbar for dispersion and intensity, default value is `true`.
% 
% `'nCol'`
% : Number of colors in the colormap, default value is 500.
% 
% `'dashed'`
% : If `true` dashed vertical lines between linear $Q$ segments will be
%   shown. Default is `false`.
% 
% `'dE'`
% : If given, a Gaussian will be convoluted with the binned data to simulate finite
%   energy resolution. Only works if `mode=3`. If zero, no convolution
%   performed. Default value is 0.
% 
% `'fontSize'`
% : Font size in pt for the labels on the plot, default value is 14 pt.
% 
% `'colormap'`
% : Colormap for plotting, default value is stored in 
%   `swpref.getpref('colormap')`. For single plot and for multiple plot it
%   will be a continuous scale from white to different color. This is the
%   `'auto'` mode. Also colormap can be given directly using standard
%   colormaps as function handles, e.g. `@jet`. To overplot multiple
%   spectra `colormap` option will be a matrix, with dimensions [3 nConv],
%   where every column defines a color for the maximum intensity. It is
%   also used for plotting dispersion curves. In case a single color all
%   dispersion curves have the same color (e.g. `[255 0 0]` for red), or as
%   many colors as dispersion curves in a matrix with dimensions of
%   $[3\times n_{mode}]$ or as a colormap function handle. In this case
%   every mode will have different color and the color is determined from
%   the index of the mode after the colormap is applied. Default value is
%   `'auto'`.
% 
% `'sortMode'`
% : Sorting the modes by energy before plotting. Default is `false`. Can 
%   improve the quality of the dispersion line plots if modes are crossing.
% 
% `'axLim'`
% : Upper limit for energy axis (for `mode` 1,2) or color axis (for `mode`
%   3), default value is `'auto'`. For color plot of multiple spectra
%   the color axis cannot be changed after the plot.
% 
% `'legend'`
% : Whether to plot legend for multiple convoluted spectras,
%   default value is `true`.
% 
% `'title'`
% : If `true` a title will be added to the figure, default value is `true`.
% 
% `'twin'`
% : Select which twins to be plotted for dispersion plots, by default the
%   spectrum corresponding to all twins will be plotted. The dimensions are
%   $[1\times n_{twinToPlot}]$.
% 
% `'lineStyle'`
% : Line style for line plots (dispersion and intensity), default value
%   `{'-' 'o-' '--'}` for plotting modes that correspond to line style of
%   $S(Q,\omega)$, $S(Q+k,\omega)$ and $S(Q-k,\omega)$ cross modes in case
%   of incommensurate magnetic systems. For commensurate systems only thte
%   first string in the cell will be considered. For example '--' gives
%   dashed lines.
% 
% `'lineWidth'`
% : Line width of line plots, default value is 0.5 pt.
% 
% `'log'`
% : If true, the 10-based logarithmic intensity will be plotted, default
%   value is `false`.
% 
% `'plotf'`
% : Function handle of the plot function for color plot. Default is
%   `surf`.
% 
% `'maxPatch'`
% : Maximum number of pixels that can be plotted using the [patch]
%   function within [sw_surf]. Using [patch] for color plot can be
%   slow on older machines, but the figure can be exported
%   afterwards as a vector graphics, using the [print] function.
%   Default value is 1000.
% 
% `'norm'`
% : If true, the convolution with a Gaussian function (in case of
%   non-zero `dE` parameter) keeps the energy integrated intensity. If
%   `false` the amplitude is kept constant. Default is determined by the
%   value stored in the input `spectra.norm`.
% 
% `'x0'`
% : Row vector with two numbers `[x0_min x0_max]`. By default the $x$ range
%   of the plot is `[0 1]` irrespective of the values of the $Q$ values. To
%   change this the lower and upper limits can be given here.
% 
% `'qlabel'`
% : Provide a list of strings for the special $Q$ points along the path in
%   the Brillouin zone, e.g. `{'X' '\Gamma' 'M' 'K' '\Gamma'}`.
% 
% `'dat'`
% : Experimental data points to plot over the calculated spectrum.
%   Can be either the name of a data file that contains the
%   experimentally fitted dispersion (needs to have the same format
%   as the input file of [spinw.fitspec] see help for details on the file
%   format), or it is a structure that contains the already imported data
%   using the [sw_readtable] function, e.g.
%
%   ```
%   T = sw_readtable('myExpData.txt','\t');
%   sw_plotspec(spectra,'dat',T);
%   ```
% 
% `'ddat'`
% : Maximum distance between any $Q$ point in the simulated spectrum
%   and an experimental data point in \\ang$^{-1}$ unit. If an
%   experimental data point is further from any $Q$ point than the given 
%   limit, it will be omitted. Default value is 0.01.
% 
% ### Output Arguments
% 
% `fHandle`
% : Handle of the plot figure.
%
% `pHandle`
% : Vector that contains the handle of the graphics objects on the figure.
% 
% ### See Also
% 
% [spinw.plot] \| [spinw.spinwave] \| [sw_surf] \| [sw_label]
%
% *[FWHM]: Full Width at Half Maximum
%

if nargin==0
    swhelp sw_plotspec
    return
end

if sw_issymspec(spectra)
    error('sw_plotspec:SymbolicInput', 'This function does not handle symbolic spectra');
end

if isfield(spectra,'norm')
    norm0 = spectra.norm;
else
    norm0 = false;
end

inpForm.fname  = {'mode' 'imag' 'aHandle' 'colorbar' 'dashed' 'norm' 'dat'     };
inpForm.defval = {4      false   gca      true       false   norm0   zeros(1,0)};
inpForm.size   = {[1 -6] [1 1]  [1 1]     [1 1]      [1 1]   [1 1]   [-9 -8]   };

inpForm.fname  = [inpForm.fname  {'dE'  'fontSize' 'colormap' 'axLim' 'ddat'}];
inpForm.defval = [inpForm.defval {0     14         'auto'     'auto'  1e-2  }];
inpForm.size   = [inpForm.size   {[1 1] [1 1]      [-1 -2]    [1 -3]  [1 1] }];

inpForm.fname  = [inpForm.fname  {'legend' 'title' 'nCol' 'twin'     'datFormat'}];
inpForm.defval = [inpForm.defval {true     true    500    zeros(1,0) 'or'       }];
inpForm.size   = [inpForm.size   {[1 1]    [1 1]   [1 1]  [1 -4]     [1 -10]    }];

inpForm.fname  = [inpForm.fname  {'lineStyle'     'lineWidth' 'sortMode'}];
inpForm.defval = [inpForm.defval {{'-' 'o-' '--'} 0.5         false     }];
inpForm.size   = [inpForm.size   {[1 -5]          [1 1]       [1 1]     }];

inpForm.fname  = [inpForm.fname  {'log' 'plotf'  'maxPatch' 'x0'  'qlabel' }];
inpForm.defval = [inpForm.defval {false @sw_surf 1000       [0 1] cell(1,0)}];
inpForm.size   = [inpForm.size   {[1 1] [1 1]    [1 1]      [1 2] [1 -7]   }];

param = sw_readparam(inpForm, varargin{:});
pref = swpref;

% plotmode string
if numel(param.mode)>1
    switch param.mode
        case 'disp'
            param.mode = 1;
        case 'int'
            param.mode = 2;
        case 'color'
            param.mode = 3;
        case {'auto' 'fancy'}
            param.mode = 4;
        otherwise
            param.mode = 4;
    end
end

% length, energy and temperature units
unitL = spectra.obj.unit.label{1};
unitE = spectra.obj.unit.label{2};
unitT = spectra.obj.unit.label{4};

% select twins for omega plot
param.twin = round(param.twin);
if isfield(spectra,'omega') && iscell(spectra.omega)
    nTwin      = numel(spectra.omega);
    if isempty(param.twin)
        param.twin = 1:nTwin;
    end
    param.twin = param.twin((param.twin<=nTwin) & (param.twin>0));
    if isempty(param.twin)
        warning('sw_plotspec:WrongInput','Number of twins is wrong, plotting all twins!');
        param.twin = 1:nTwin;
    end
    nTwin      = numel(param.twin);
else
    nTwin = 1;
    param.twin = 1;
end

if ~isfield(spectra,'omega')
    param.mode = 3;
end

if ~isfield(spectra,'swConv') && param.mode>1 && param.mode<4
    error('sw_plotspec:WrongInput',['Reference to non-existent field ''swConv'','...
        'use ''sw_egrid'' to produce the convoluted spectra before plotting!'])
end

% select twins for convoluted plots
if param.mode>1 && param.mode<4 && iscell(spectra.swConv)
    % number of convoluted spectras to plot
    nTwinS      = size(spectra.swConv,2);
    param.twinS = param.twin((param.twin<=nTwinS) & (param.twin>0));
    if isempty(param.twinS) && (nTwinS>1)
        warning('sw_plotspec:WrongInput','Number of twins is wrong, plotting all twins!');
        param.twinS = 1:nTwinS;
    end
    nTwinS      = numel(param.twinS);
else
    nTwinS      = 1;
    param.twinS = 1;
end

% Determine powder mode
powmode = false;
if numel(spectra.hklA)==length(spectra.hklA)
    powmode = true;
end

if param.mode == 4
    % PLOT EASY PEASY
    
    fHandle = [];
    pHandle = [];
    pColor = isfield(spectra,'swConv');
    
    if pColor
        if param.dE == 0
            Eres = (spectra.Evect(end) - spectra.Evect(1))/50;
        else
            Eres = param.dE;
        end
        
        [fHandle, pHandle] = sw_plotspec(spectra,'mode',3,'dE',Eres,...
            'dashed',true,'colorbar',false,'axLim',param.axLim,...
            'lineStyle',param.lineStyle,'maxPatch',...
            param.maxPatch,'qLabel',param.qlabel,'dat',param.dat,...
            'ddat',param.ddat,'datFormat',param.datFormat,...
            'legend',param.legend);
    end
    if ~powmode
        hold on
        if pColor
            cMap0 = [0 0 0];
        else
            cMap0 = 'auto';
        end
        
        if iscell(spectra.omega)
            omegaTemp = cell2mat(spectra.omega);
            Emax = max(real(omegaTemp(:)));
            clear('omegaTemp');
        else
            Emax = max(real(spectra.omega(:)));
        end

        [fHandle, pHandle] = sw_plotspec(spectra,'mode','disp','colorbar',~pColor,...
            'dashed',false,'title',~pColor,'legend',~pColor,'imag',~pColor,...
            'lineStyle',param.lineStyle,'colormap',cMap0,'axLim',[0 1.1*Emax],...
            'qLabel',param.qlabel);
    end
    
    if nargout >0
        fHandle0 = fHandle;
    end
    if nargout>1
        pHandle0 = pHandle;
    end
    return
end

% Label of the x-axis
if powmode
    % powder mode
    xLabel  = ['Momentum transfer (' unitL '^-1)'];
    xAxis   = spectra.hklA;
else
    [xLabel, xAxis] = sw_label(spectra.hkl,spectra.hklA,spectra.obj.unit.label{1});
    if ~isempty(param.qlabel) && iscell(xLabel)
        if numel(param.qlabel)~=(numel(xLabel)-1)
            error('sw_plotspec:WrongInput','The number of q labels is wrong!')
        end
        % change labels
        xLabel(1:(end-1)) = param.qlabel;
        xLabel0 = 'Momentum';
    else
        xLabel0 = 'Momentum (r.l.u.)';
    end
    
    % shift the x-axis if requested
    xAxis = xAxis*diff(param.x0)+param.x0(1);
end

if isfield(spectra,'Evect')
    % create center bin
    Evect = spectra.Evect;
    yAxis = (Evect(2:end)+Evect(1:(end-1)))/2;
else
    yAxis = [min(spectra.omega(:)) max(spectra.omega(:))];
end

yLabel = ['Energy transfer (' unitE ')'];

if isa(param.aHandle,'matlab.graphics.axis.Axes') || ishandle(param.aHandle) 

    if any(strfind(get(get(param.aHandle,'Parent'),'Tag'),'sw_crystal'))
        % don't plot into the crystal structure window
        fHandle = figure('color','w');
        param.aHandle = gca;
    else
        fHandle = get(gca,'Parent');
    end
    axes(param.aHandle);
else
    fHandle = sw_getfighandle('sw_spectra');
    if isempty(fHandle)
        fHandle = figure('color','w');
    end
end
% set Tag to find window later easily
set(fHandle,'Tag','sw_spectra');

setappdata(fHandle,'param',param);
setappdata(fHandle,'spectra',spectra);


% Plotting styles for commensurate/incommensurate structures.
if iscell(param.lineStyle)
    if numel(param.lineStyle) ~= 3
        param.lineStyle = repmat(param.lineStyle(1),[1 3]);
    end
else
    param.lineStyle = repmat({param.lineStyle},[1 3]);
end


if ~powmode
    % sort the convoluted intensities into cell array.
    if param.mode>1
        if ~iscell(spectra.component)
            swInt     = {spectra.swInt};
            swConv    = {spectra.swConv};
            component = {spectra.component};
        else
            swInt     = spectra.swInt(:,param.twinS);
            swConv    = spectra.swConv(:,param.twinS);
            component = spectra.component;
        end
        % number of different convoluted cross sections
        nConv = numel(component);
        % number of convoluted plots
        nPlot  = nTwinS * nConv;
    end
    
    nMagExt = spectra.obj.nmagext;
    % package all fields into cells for easy looping over twins
    if isfield(spectra,'omega')
        if ~iscell(spectra.omega)
            omega = {spectra.omega};
        else
            omega = spectra.omega(:,param.twin);
        end
        nMode   = size(omega{1},1);
    else
        omega = {0};
        nMode = 2*nMagExt;
    end
    
    
    if param.mode<3
        % Defines colors for plotting modes.
        %colors  = flipud(fireprint(nMode+2));
        if isa(param.colormap,'function_handle')
            colors = flipud(param.colormap(nMode+2));
        else
            if strcmpi(param.colormap,'auto')
                colors = flipud(cm_fireprint(nMode+2));
            else
                if numel(param.colormap) == 3
                    param.colormap = param.colormap(:);
                    colors = repmat(param.colormap',nMode+2,1)/255;
                elseif (size(param.colormap,1) == nMode) && (size(param.colormap,2)==3)
                    colors = [0 0 0; param.colormap; 0 0 0]/255;
                elseif (size(param.colormap,2) == nMode) && (size(param.colormap,1)==3)
                    colors = [0 0 0; param.colormap'; 0 0 0]/255;
                else
                    error('sw_plotspec:ColormapError','The dimensions of the colormap should be [3 nMode=%d]',nMode);
                end
            end
        end
        colors  = colors(2:(end-1),:);
    end
    
    modeList = nMode/(2*nMagExt);
    if modeList == 1
        if param.imag
            lLabel = {'Real' 'Imaginary'};
        else
            lLabel = {'Real'};
        end
    elseif modeList == 3
        if param.imag
            lLabel = {'Q+k_m' 'Q' 'Q-k_m' 'Imaginary'};
        else
            lLabel = {'Q+k_m' 'Q' 'Q-k_m'};
        end
    else
        error('sw_plotspec:NumberOfModes','Wrong number of spin wave modes!');
    end
else
    nPlot     = 1;
    swConv    = {spectra.swConv};
    nConv     = 1;
    component = {spectra.component};
    param.legend = false;
end

if powmode && (param.mode~=3)
    warning('sw_plotspec:PowMode','Powder spectra, only convoluted spectra can be plotted!');
    param.mode = 3;
end

hPlot = gobjects(0);
hold on

switch param.mode
    case 1
        % Line plot of dispersion
        axis0 = [xAxis(1) xAxis(end) yAxis(1) yAxis(end)];
        titleStr0 = 'Spin wave dispersion: \omega(Q)';
        % loop over the twins
        for tt = 1:nTwin
            plotr = (real(omega{1,tt}));
            ploti = (imag(omega{1,tt}));
            if param.sortMode
                plotr = sort(plotr,1);
                ploti = sort(ploti,1);
            end
            % loop over all spin wave modes
            for ii = 1:nMode
                incIdx = ceil(ii/2/nMagExt);
                hPlot(end+1)    = plot3(xAxis,plotr(ii,:),xAxis*0+1e5,param.lineStyle{incIdx},...
                    'Color', colors(ii,:),'LineWidth',param.lineWidth); %#ok<*AGROW>
                hLegend(incIdx) = hPlot(end);
                if param.imag
                    hPlot(end+1)        = plot(xAxis,ploti(ii,:),'ro-');
                    hLegend(modeList+1) = hPlot(end);
                end
            end
        end
        
    case 2
        % Line plot of cross sections but only the first cell array element
        axis0 = [xAxis(1) xAxis(end) 0 1];
        
        yLabel = 'Intensity (arb. u.)';
        if param.log
            yLabel = ['log ' yLabel];
        end
        
        titleStr0 = 'Intensity of the spin-spin correlation function: ';
        % loop over the twins
        for tt = 1:nTwinS
            for jj = 1:nConv
                if param.imag
                    plotr = abs(imag(swInt{jj,tt}));
                else
                    plotr = abs(real(swInt{jj,tt}));
                end
                
                if param.log
                    plotr = log10(plotr);
                end
                
                for ii = 1:nMode
                    hPlot(end+1) = plot3(xAxis,plotr(ii,:),xAxis*0+1e5,param.lineStyle{mod(jj-1,3)+1},...
                        'Color', colors(ii,:),'LineWidth',param.lineWidth); %#ok<*AGROW>
                    if ii == nMode
                        hLegend(jj) = hPlot(end);
                    end
                end
            end
            if param.imag
                
                for jj = 1:nConv
                    ploti = abs(imag(swInt{jj,tt}));
                    for ii = 1:nMode
                        hPlot(end+1) = plot3(xAxis,ploti(ii,:),xAxis*0+2e5,'ro-');
                        hLegend(nConv+1) = hPlot(end);
                    end
                end
                
            end
        end
end

if param.mode < 3
    if strcmpi(param.axLim,'auto')
        autAxis = axis;
        axis([sort(axis0(1:2)) autAxis(3:4)]);
    else
        if numel(param.axLim) == 1
            param.axLim = [0 param.axLim];
        end
        axis([sort(axis0(1:2)) param.axLim]);
    end
    box on
end

if param.mode == 1
    if param.legend
        legend(hLegend,lLabel{:});
    end
    
    if param.colorbar
        cHandle = colorbar;
        set(get(cHandle,'ylabel'),'String', 'Index of dispersion line (each color different mode)');
        colormap(colors);
        caxis([0.5 nMode+0.5]);
        set(cHandle,'YTick',1:nMode);
    end
end

% current axis
hAxis = gca;

if param.mode == 3
    
    % filter out imaginary, inf and NaN values
    mask = cell(1,nPlot);
    for ii = 1:nPlot
        if param.imag
            swConv{ii} = imag(swConv{ii});
        else
            swConv{ii} = real(swConv{ii});
        end
        maskT = swConv{ii};
        maskT(~isnan(maskT)) = 1;
        mask{ii} = maskT;
        
        swConv{ii}(isnan(swConv{ii})) = 0;
        swConv{ii}(isinf(swConv{ii})) = 0;
    end
    
    % make colormap from white to given color
    % by defining a cell of functions, with dimension of [1 nPlot]
    if ~isa(param.colormap,'function_handle')
        if strcmpi(param.colormap,'auto')
            % for 'auto' mode an equally space hue values are created for
            % use with multiple sectra plot
            if nPlot>1
                param.colormap = hsv2rgb([(1:nPlot)'/nPlot ones(nPlot,2)])'*255;
            else
                param.colormap = {pref.colormap};
            end
        end
        if ~iscell(param.colormap)
            if (size(param.colormap,1) ~= 3) || (size(param.colormap,2)<nPlot)
                error('sw_plotspec:ColormapError','The dimensions of the colormap should be [3 nPlot]');
            end
            tHandle = cell(1,nPlot);
            for ii = 1:nPlot
                tHandle{ii} = @(numSteps)makecolormap(param.colormap(:,ii)'/255,[1 1 1],numSteps);
            end
            param.colormap = tHandle;
        end
    else
        param.colormap = repmat({param.colormap},[1 nPlot]);
    end
    
    % Gaussian energy convolution kern
    if param.dE>0
        sG = param.dE/2.35482;
        %x0 = spectra.Evect;
        x0 = yAxis;
        dx = (x0(end)-x0(1))/(length(x0)-1);
        nG = ceil(3*sG/dx);
        xG = (-nG*dx):dx:(nG*dx);
        % Gaussian normalised intensity
        fG = exp(-(xG/sG).^2/2);
        if param.norm
            fG = fG/sum(fG);
        else
            fG = fG/max(fG(:));
        end
    else
        fG = 1;
    end
    
    % c axis limit
    if strcmpi(param.axLim,'auto')
        % determine maximum intensity automatically
        zi = reshape(cell2mat(swConv),1,[]);
        
        zi(isnan(zi)) = [];
        posLim = sort([zi(zi>0) 0],'descend');
        negLim = sort([zi(zi<0) 0],'ascend');
        axLim  = [negLim(ceil(end*8e-2)) posLim(ceil(end*8e-2))];
        axMagn = 10.^(floor(log10(abs(axLim))));
        axLim  = ceil(axLim./axMagn).*axMagn;
        axLim(isnan(axLim)) = 0;
        cMaxMax = max(abs(zi));
        if axLim(2)-axLim(1) == 0
            axLim = [0 1];
        end
        
        if axLim(1) < -1e-8
            axLim = [-max(abs(axLim)) max(abs(axLim))];
        end
        %         zi(zi<=0) = [];
        %         maxsort  = sort(zi,'descend');
        %         maxsortS = maxsort(ceil(end*8e-3));
        %         magni    = 10^(floor(log10(maxsortS)));
        %         cMax     = ceil(maxsortS/magni)*magni;
        %         cMaxMax  = maxsort(1);
        %
        %         if cMax <= 0
        %             cMax = 1;
        %         end
    else
        axLim   = param.axLim;
        cMaxMax = max(abs(axLim));
        
        if numel(axLim) == 1
            axLim = [0 axLim];
        end
    end
    
    % convolute spectra with Gaussian to simulate finite energy resolution
    for ii = 1:nPlot
        iTemp = swConv{ii}';
        iTemp = conv2(iTemp,fG,'same');
        swConv{ii} = iTemp';
        % take log if param.log is true
        if param.log
            swConv{ii} = log10(swConv{ii}+1e-10);
        end
    end
    
    
    if nPlot == 1
        % single spectra
        imageDisp = (swConv{1}.*mask{ii})';
        
        % Use surf to hide the NaN numbers
        [X, Y] = meshgrid(xAxis,yAxis);
        cMap = flipud(param.colormap{1}(param.nCol));
        
        if cMaxMax <1e-6
            %hSurf = param.plotf(X,Y,imageDisp'*0);
            axLim = [0 1];
            param.plotf(X,Y,imageDisp'*0,axLim,cMap,param.maxPatch);
            
        else
            % hSurf = param.plotf(X,Y,imageDisp');
            param.plotf (X,Y,imageDisp',axLim,cMap,param.maxPatch);
        end
        %view(2);
        %if all(ishandle(hSurf))
        %    set(hSurf,'EdgeAlpha',0);
        %end
        
        
    else
        
        % multiple spectra
        vMat = zeros([size(swConv{1}),nPlot]);
        for ii = 1:nPlot
            vMat(:,:,ii) = swConv{ii};%.*mask{ii};
        end
        cMat = sw_multicolor(vMat, param.colormap, axLim, param.nCol,true);
        % plot image piece-by-pice for the different Q directions
        if iscell(xLabel)
            xCut  = xLabel{end};
            nCut  = numel(xCut);
            hPlot = gobjects(1,nCut);
            for ii = 2:nCut
                selIdx = xAxis>=xCut(ii-1) & xAxis<=xCut(ii);
                hPlot(ii-1) = image(xAxis(selIdx),yAxis,cMat(:,selIdx,:));
            end
            
        else
            hPlot = image(xAxis,yAxis,cMat);
        end
        set(gca,'YDir','normal');
    end
    
    caxis(axLim);
    
    if param.colorbar && (nPlot == 1)
        cHandle = colorbar;
        if spectra.norm
            if spectra.obj.unit.nformula > 0
                cLabelU = ['(mbarn/' unitE '/f.u.)'];
            else
                cLabelU = ['(mbarn/' unitE '/cell)'];
            end
        else
            cLabelU = '(arb. u.)';
        end
        if param.log
            cLabel = ['log Intensity ' cLabelU];
        else
            cLabel = ['Intensity ' cLabelU];
        end
        
        set(get(cHandle,'ylabel'),'String',cLabel);
    end
    for ii = 1:nPlot
        lColor      = param.colormap{ii}(2);
        hLegend(ii) = plot(-1e5,-1e5,'Color',lColor(1,:));
    end
    
    axis0 = [xAxis(1) xAxis(end) yAxis(1) yAxis(end)];
    axis([sort(axis0(1:2)) axis0(3:4)])
    
    grid off
    set(gca,'Layer','top')
    titleStr0 = 'Convoluted ';
    if powmode
        titleStr0 = [titleStr0 'powder '];
    end
    titleStr0 = [titleStr0 'spectra: '];
    box on
    
    % overplot data
    if ~isempty(param.dat)
        % read table of experimental data
        if ischar(param.dat) || iscell(param.dat)
            if ~iscell(param.dat)
                param.dat = {param.dat};
            end
            T = sw_readtable(param.dat{:});
        else
            T = param.dat;
        end
        
        % collect experimental data into matrices
        Qexp = [[T(:).QH];[T(:).QK];[T(:).QL]];
        nQ = size(Qexp,2);
        % number of modes
        nMode = sum(cellfun(@(C)numel(C)>1 && strcmp(C(1),'I') ,fieldnames(T)));
        
        iName = strsplit(strtrim(sprintf('I%d ',1:nMode)),' ');
        eName = strsplit(strtrim(sprintf('EN%d ',1:nMode)),' ');
        sName = strsplit(strtrim(sprintf('s%d ',1:nMode)),' ');
        
        dat.I = zeros(nMode,nQ);
        dat.E = zeros(nMode,nQ);
        dat.s = zeros(nMode,nQ);
        
        for ii = 1:nMode
            dat.I(ii,:) = [T.(iName{ii})];
            dat.E(ii,:) = [T.(eName{ii})];
            dat.s(ii,:) = [T.(sName{ii})];
        end
        
        dat.s(dat.I==0) = nan;
        dat.E(dat.I==0) = nan;
        % reciprocal lattice
        RL   = spectra.obj.rl;
        
        % distance of experimental data points from plotted data points
        D = sqrt(sum(bsxfun(@minus,permute(Qexp'*RL,[1 3 2]),permute(spectra.hkl'*RL,[3 1 2])).^2,3));
        % idx stores the index of the point
        [sel,idxD] = min(D,[],2);
        % experimental data points that will appear on the plot
        sel = sel < param.ddat;
        idxD = idxD(sel);
        
        % add new axis
        hAxis(2) = axes('Position',hAxis.Position,'Color','none');
        linkaxes(hAxis,'xy');
        hold on
        if ~iscell(param.datFormat)
            param.datFormat = {param.datFormat};
        end
        
        for jj = 1:nMode
            errorbar(xAxis(idxD),dat.E(jj,sel),dat.s(jj,sel),param.datFormat{:})
        end
        axes(hAxis(1));
    end
    
end

ylabel(yLabel);

if param.mode > 1
    titleStr = cell(1,nConv);
    for ii = 1:nConv
        titleStr{ii} = sw_titlestr(component{ii});
        if param.imag
            titleStr{ii} = ['Im ' titleStr{ii}];
        else
            titleStr{ii} = ['Re ' titleStr{ii}];
        end
        
    end
    
    for ii = 1:nConv-1
        titleStr0 = [titleStr0 titleStr{ii} ', '];
    end
    titleStr0 = [titleStr0 titleStr{end}];
    
    if param.legend
        % put the names of the twins into the plot
        titleStr = repmat(titleStr,[1 nTwinS]);
        if nTwinS > 1
            for ii = 1:nPlot
                titleStr{ii} = ['Twin #' num2str(ceil(ii/nConv)) ' ' titleStr{ii}];
            end
        end
        % put label for imaginary lineplot
        if param.imag && param.mode==2
            titleStr = [titleStr(1:nConv) lLabel(2)];
        end
        warn_state = warning;
        warning('off','MATLAB:legend:IgnoringExtraEntries')
        set(legend(hLegend,titleStr{:}),'FontSize',param.fontSize);
        warning(warn_state);
    end
end

if param.title
    if isfield(spectra,'T')
        title([titleStr0 sprintf([', T = %5.1f ' unitT],spectra.T)],'FontSize',param.fontSize);
    else
        title(titleStr0,'FontSize',param.fontSize);
    end
end
drawnow;
hold off

if iscell(xLabel)
    xTickLoc = xLabel{end};
    set(gca,'XTick',xTickLoc)
    set(gca,'XTickLabel',xLabel(1:end-1));
    if param.dashed
        autAxis = axis;
        for jj=2:length(xLabel)-1
            line([1 1]*xTickLoc(jj),autAxis(3:4),'LineStyle','--','color',[0 0 0]);
        end
    end
    xlabel(xLabel0);
else
    xlabel(xLabel);
end

if nargout == 1
    fHandle0 = fHandle;
elseif nargout == 2
    fHandle0 = fHandle;
    pHandle0 = hPlot;
end

if numel(hAxis)>1
    legend off
    axes(hAxis(2));
    hAxis(2).Visible = 'off';
    hAxis(2).Position = hAxis(1).Position;
end

end

function titleStr = sw_titlestr(component)
% creates plot title string from param.component string

titleStr = [component '}(\omega,Q)'];
titleStr = strrep(titleStr,'perp','\perp ');
titleStr = strrep(titleStr,'S','S^{');
titleStr = strrep(titleStr,'P','P^{');
titleStr = strrep(titleStr,'M','M^{');
titleStr = strrep(titleStr,'+','}(\omega,Q) +');
titleStr = strrep(titleStr,'-','}(\omega,Q) -');

end

function [xLabel, xAxis] = sw_label(hkl,hklA,lUnit)
% returns axis labels for spectrum plot
% 
% ### Syntax
% 
% `[xlabel, xaxis] = sw_label(hkl,hkla)`
% 
% ### Description
% 
% `[xlabel, xaxis] = sw_label(hkl,hkla)` returns the label for the x-axis
% and x-coordinates for a 
% 
% ### Input Arguments
% 
% `hkl`
% : Momentum transfer values in r.l.u., dimensions are [3 nQ].
% 
% `hklA`
% : Momentum transfer values in \\ang$^{-1}$, dimensions are [3 nQ].
% 
% `lUnit`
% : Length unit, given in a string.
% 
% ### Output Arguments
% 
% It returns the label and axis vector for the x-axis for momentum transfer
% scans linear in reciprocal space.
% 
% ### See Also
% 
% [sw_plotspec]
%

if nargin == 0
    help sw_label
    return
end

if nargin<3
    lUnit = symbol('a');
end

hkl  = hkl';
hklA = hklA';

nQ  = size(hkl,1);

% distance between start and end points
dk0 = hkl(1,:) - hkl(end,:);

% determine whether it is line scan
if abs(dk0(1))>1e-5
    dk    = dk0/dk0(1);
    xAxis = hkl(:,1);
elseif abs(dk0(2))>1e-5
    dk    = dk0/dk0(2);
    xAxis = hkl(:,2);
elseif abs(dk0(3))>1e-5
    dk    = dk0/dk0(3);
    xAxis = hkl(:,3);
end

% parse curve into straight lines
qStep  = hkl(2:end,:)-hkl(1:(end-1),:);
qStep  = bsxfun(@rdivide,qStep,sqrt(sum(qStep.^2,2)));
qCurve = sum(qStep(2:end,:).*qStep(1:end-1,:),2);
qIdx   = find(qCurve<0.97)+1;

if numel(qIdx) == 0
    linescan = 0;
elseif numel(qIdx)/nQ<0.2
    % line scan with straight pieces
    linescan = 1;
else
    % curved scan
    linescan = 2;
end

switch linescan
    case 0
        % single linear scan in r.l.u.
        k0 = hkl(1,:);
        
        changeX = false;
        inA = sqrt(sum(((hklA(1,:)' - hklA(end,:)').^2)))/abs(xAxis(end)-xAxis(1));
        xiLabel = cell(1,3);
        for ii = 1:3
            if abs(k0(ii)) > 1e-5
                if abs(dk(ii)) > 1e-3
                    if abs(dk(ii)-1)<1e-3
                        xiLabel{ii} = sprintf('%.4g+\\xi',k0(ii));
                        if ~changeX
                            xAxis = xAxis - k0(ii);
                            changeX = true;
                        end
                    elseif abs(dk(ii)+1)<1e-3
                        xiLabel{ii} = sprintf('%.4g-\\xi',k0(ii));
                        if ~changeX
                            xAxis = xAxis + k0(ii);
                            changeX = true;
                        end
                        
                    else
                        xiLabel{ii} = sprintf('%.4g%+.4g\\xi',k0(ii),dk(ii));
                        if ~changeX
                            xAxis = xAxis - k0(ii);
                            changeX = true;
                        end
                        
                    end
                    
                else
                    xiLabel{ii} = sprintf('%.4g',k0(ii));
                end
            else
                if abs(dk(ii)) > 1e-3
                    if abs(dk(ii)-1)<1e-3
                        xiLabel{ii} = '\xi';
                    elseif abs(dk(ii)+1)<1e-3
                        xiLabel{ii} = '-\xi';
                    else
                        xiLabel{ii} = sprintf('%.4g\\xi',dk(ii));
                    end
                else
                    xiLabel{ii} = '0';
                end
            end
        end
        if size(hkl,1) == 1
            xLabel = ['(' xiLabel{1} ',' xiLabel{2} ',' xiLabel{3} ')'];
        else
            xLabel = ['(' xiLabel{1} ',' xiLabel{2} ',' xiLabel{3} ') in ' sprintf('%.5g ',inA) lUnit '^{-1}'];
        end
        
    case 1
        % use inverse Angstrom for the x-axis scaling
        qIdx = [1;qIdx;nQ];
        xAxis = 0;
        for ii = 2:length(qIdx)
            hklAdist = sqrt(sum((hklA(qIdx(ii),:)-hklA(qIdx(ii-1),:)).^2));
            qAdd = linspace(0,hklAdist,qIdx(ii)-qIdx(ii-1)+1);
            qAdd = qAdd(2:end);
            xAxis = [xAxis qAdd+xAxis(end)];
        end
        % create labels for line pieces
        xLabel = cell(1,length(qIdx));
        for ii = 1:length(qIdx)
            hkl(abs(hkl)<1e-3) = 0;
            xLabel{ii} = ['(' num2str(hkl(qIdx(ii),1)) ',' num2str(hkl(qIdx(ii),2)) ',' num2str(hkl(qIdx(ii),3)) ')'];
        end
        xLabel{end+1} = xAxis(qIdx);
        1;
    case 2
        xLabel = 'Momentum transfer';
        xAxis  = linspace(0,1,nQ);
end

end
