function spectra = sw_egrid(spectra, varargin)
% calculates energy bins of a spectrum 
% 
% ### Syntax
% 
% `spectra = sw_egrid(spectra,Name,Value)`
% 
% ### Description
% 
% `spectra = sw_egrid(spectra,Name,Value)` takes a calculated spectrum that
% contains $S^{\alpha\beta}(Q,\omega)$ and converts it into an intensity
% map `I(i,j)` via binning the energy values and selecting a given
% component of the $9\times 9$ spin-spin correlation matrix. For example by
% default (setting the `component` parameter to `'Sperp'`) it selects the
% neutron scattering cross section via calculating the following quantity:
%
%   $S_\perp(Q,\omega)=\sum_{\alpha\beta}(1-\hat{q}^\alpha\hat{q}^\beta)\cdot S^{\alpha\beta}(Q,\omega)$
%  
% 
% ### Examples
% 
% The line will create an energy bin, with steps of 0.1 and bins the
% spin-spin correlation function. Two different matrices will be
% calculated, first using the sum of the $S^{xx}$ and $S^{yy}$ components, second
% will contain the $S^{zz}$ component of the correlation function.
%
% ```
% >>tri = sw_model('triAF',1)
% >>spectra = tri.spinwave({[0 0 0] [1 1 0] 501})
% >>E = linspace(0,5,501)
% >>spectra = sw_egrid(spectra,'component',{'Sxx+Syy' 'Szz'},'Evect',E)
% >>figure
% >>sw_plotspec(spectra,'mode','color','axLim',[0 0.5],'dE',0.2)
% >>snapnow
% ```
%
% ### Input Arguments
% 
% `spectra`
% : Input structure, contains spin-spin correlation functions. Supported
%   inputs are produced by [spinw.spinwave], [spinw.powspec] and
%   [spinw.scga].
% 
% ### Name-Value Pair Arguments
% 
% `'component'`
% : A string that Selects a correlation function component that will be
%   binned. The possible values are:
%   * `'Sperp'` bins the magnetic neutron scattering intensity
%     (the $\langle S_\perp S_\perp\rangle$ expectation value). Default.
%   * `'Sab'`   bins the selected components of the spin-spin
%               correlation function. Letter `a` and `b` can be `x`,
%               `y` or `z`. For example: `'Sxx'` will convolute the
%               $S^{xx}(Q,\omega)$ component of the correlation function with the
%               dispersion. Here the $xyz$ is the standard coordinate system.
%   *`'Mab'`    bins the selected components of the spin-spin
%               correlation function in the Blume-Maleev coordinate system.
%               Letter `a` and `b` can be `x`, `y` or `z`. For example:
%               `'Mxx'` will convolute the `xx` component of the
%               correlation function with the dispersion.
%   * `'Pab'`   bins the selected component of the polarisation
%               matrix. Letter `a` and `b` can be `x`, `y` or `z`. For
%               example: `'Pyy'` will convolute the `yy` component of
%               the polarisation matrix with the dispersion. The
%               coordinates used are in the Blume-Maleev coordinate
%               system, see below.
%   * `'Pa'`    bins the intensity of the calculated polarised
%               neutron scattering, with inciden polarisation of
%               `Pa` where letter `a` can be `x`, `y` or `z`. For example:
%               `'Py'` will convolute the scattering intensity
%               simulated for incident polarisation $P_i\|y$. The
%               used coordinates are in the Blume-Maleev coordinate
%               system.
%   * `'fName'` where `fName` is one of the field names of the input
%               structure spectra. This field should contain a
%               matrix with dimensions of $[n_{mode}\times n_{hkl}]$.
%
%   Any linear combination of the above are allowed, for example:
%   `'Sxx+2*Syy'` will bin the linear combination of the `xx` component of
%   the spin-spin correlation function with the `yy` component.
%   Several cross section can be convoluted and stored
%   independently, if component is a cell array containing strings
%   each containing any linear combination of cross sections as
%   above, the cell array needs to have size $[1\times n_{cell}]$, for
%   example `{'Sxx' 'Syy' 'Szz'}`.
% 
% `'Evect'`
% : Row vector that defines the center/edge of the energy bins of the
%   calculated output, number of elements is $n_E$. The energy units
%   are defined by the [spinw.unit] property. Default
%   value is an edge bin: `linspace(0,1.1*maxOmega,501)`.
% 
% `'binType'`
% : String, determines the type of bin give, possible options:
%   * `'cbin'`      Center bin, the center of each energy bin is given.
%   * `'ebin'`      Edge bin, the edges of each bin is given.
%   Default value is `'ebin'`.
% 
% `'T'`
% : Temperature, used to calculate the Bose factor in units
%   depending on the Boltzmann constant stored in [spinw.unit]. Default
%   temperature is taken from `obj.single_ion.T`. The Bose factor is
%   included in `swConv` field of the output.
% 
% `'sumtwin'`
% : If true, the spectra of the different twins will be summed
%   together weighted with the normalized volume fractions, see
%   [spinw.twin]. Default value is true.
% 
% `'modeIdx'`
% : Select certain spin wave modes from the $2*n_{magatom}$ number of
%   modes to include in the output. Default value is `1:2*nMagAtom` to
%   include all modes.
% 
% `'epsilon'`
% : DEPRECATED (previously the error limit, used to determine whether a 
%   given energy bin is uniform or not. Default value is $10^{-5}$). This
%   parameter is no longer relevant and is ignored.
% 
% `'autoEmin'`
% : Due to the finite numerical precision, the spin wave energies
%   can contain small imaginary values. These can ruin the
%   convoluted spectrum at low energies. To improve the spectrum,
%   the lowest energy bin should start above the imaginary part of
%   the spin wave energy. If `'autoEmin'` is set to `true`, it
%   calculates the bottom of the first energy bin automatically and
%   overwrites the given value. Only works if the input energy bin
%   starts with zero. Default value is `false`.
% 
% `'imagChk'`
% : Checks whether the imaginary part of the spin wave dispersion is
%   smaller than the energy bin size. Default value is true.
% 
% {{note The Blume-Maleev coordinate system is a cartesian coordinate
% system with $x_{BM}$, $y_{BM}$ and $z_{BM}$ basis vectors defined as:
% <br> $x_{BM}$    parallel to the momentum transfer $Q$,
% <br> $y_{BM}$    perpendicular to $x_{BM}$ in the scattering plane,
% <br> $z_{BM}$    perpendicular to the scattering plane.
% }}
% 
% ### Output Arguments
% 
% `spectra` same as the input `spectra` plus additions fields:
% 
% `swConv`
% : Stores the selected cross section binned in energy in a matrix with
%   dimensions of $[n_E\times n_{hkl}]$. Includes the Bose factor.
% 
% `swInt`
% : Stores the selected cross sections for every mode in a matrix with
%   dimensions of $[n_{mode}\times n_{hkl}]$.
% 
% `T`
% : Input temperature.
% 
% `component`
% : Cell that contains the input component selector strings.
% 
% `Evect`
% : Input energy bin vector, defines the energy bin **edge** positions
%   (converted from the given bin centers if necessary).
%
% `zeroEnergyTol`
% : Eigenvalues with magnitude of the real component less than zeroEnergyTol
%   will be not be included in the structure factor binning
% 
% `param`
% : All the input parameters.
%
% If `'component'` parameter is a cell array or the spectra of multiple
% twins are convoluted separately, swConv and swInt will be a cell that
% packages the matrices corresponding to each component/twin. The
% dimensions of the cell are $[n_{conv}\times n_{twin}]$.
% 
% ### See Also
% 
% [spinw.spinwave] \| [sw_neutron]
%

if nargin == 0
    swhelp sw_egrid
    return
end

if sw_issymspec(spectra)
    error('sw_egrid:SymbolicInput', 'This function does not handle symbolic spectra');
end

if isfield(spectra,'obj')
    T0 = spectra.obj.single_ion.T;
else
    T0 = 0;
end

% use Evect if defined in spectra
if isfield(spectra,'Evect')
    E0 = spectra.Evect;
else
    E0 = [];
end

inpForm.fname  = {'Evect' 'T'   'component' 'sumtwin' 'modeIdx'  'epsilon'};
inpForm.defval = {E0      T0    'Sperp'     true      zeros(1,0) NaN     };
inpForm.size   = {[1 -1]  [1 1] [1 -2]      [1 1]     [1 -4]     [1 1]    };
inpForm.soft   = {true    false false       false     false      false    };

inpForm.fname  = [inpForm.fname  {'autoEmin' 'imagChk' 'binType' 'zeroEnergyTol'}];
inpForm.defval = [inpForm.defval {false       true       'ebin'  5e-4}];
inpForm.size   = [inpForm.size   {[1 1]      [1 1]      [1 -5]   [1 1]}];
inpForm.soft   = [inpForm.soft   {false      false      false    false}];

inpForm.fname  = [inpForm.fname  {'maxDSF'}];
inpForm.defval = [inpForm.defval {1e6}];
inpForm.size   = [inpForm.size   {[1 1]}];
inpForm.soft   = [inpForm.soft   {false}];

param = sw_readparam(inpForm, varargin{:});

if ~isnan(param.epsilon)
    % non default value of epsilon
    warning('sw_egrid:DeprecationWarning', ...
            'epsilon is deprecated - it is not relevant and will be ignored.');
end

switch param.binType
    case 'ebin'
        eBin = true;
        param.binType = 1;
    case 'cbin'
        eBin = false;
        param.binType = 2;
    otherwise
        error('sw_egrid:WrongInput','Wrong energy bin type!')
end

if isempty(param.Evect)
    if isfield(spectra,'omega')
        if iscell(spectra.omega)
            omegaTemp = cell2mat(spectra.omega);
            Emax = max(real(omegaTemp(:)));
            clear('omegaTemp');
        else
            Emax = max(real(spectra.omega(:)));
        end
        if Emax == 0
            Emax = 1;
        end
        
        param.Evect = linspace(0,1.1*Emax,501);
        eBin = true;
    else
        param.Evect = [-1 1];
        eBin = true;
    end
end

% parse the component string
if iscell(param.component)
    nConv = numel(param.component);
    parsed = cell(1,nConv);
    for ii = 1:numel(param.component)
        parsed{ii} = sw_parstr(param.component{ii});
    end
elseif isstruct(param.component)
    nConv  = 1;
    parsed = {param.component};
    param.component = {parsed{1}.string};
else
    nConv = 1;
    parsed = {sw_parstr(param.component)};
    param.component = {param.component};
end

% get all the used type of correlation functions
parType = [];
for ii = 1:numel(parsed)
    for jj = 1:numel(parsed{ii}.type)
        parType = [parType parsed{ii}.type{jj}(1)]; %#ok<AGROW>
    end
end
parType = ismember(1:5,parType);

if parType(1)
    if ~isfield(spectra,'Sperp') || isempty(spectra.Sperp)
        % auto produce neutron scattering cross section
        %fprintf0(fid,'Neutron scattering cross section is calculated.\n');
        spectra = sw_neutron(spectra);
    end
    Sperp = spectra.Sperp;
else
    Sperp = [];
end

if parType(3)
    if isfield(spectra,'Mab') && ~isempty(spectra.Mab)
        Mab = spectra.Mab;
    else
        error('sw_egrid:WrongInput',['Reference to non-existent field ''Mab'','...
            ' use ''sw_neutron(''pol'',true)'' to produce the polarised neutron scattering cross sections,'...
            ' before binning with energy transfer!'])
    end
else
    Mab = [];
end

if parType(4)
    if isfield(spectra,'Pab') && ~isempty(spectra.Pab)
        Pab = spectra.Pab;
    else
        error('sw_egrid:WrongInput',['Reference to non-existent field ''Pab'','...
            ' use ''sw_neutron(''pol'',true)'' to produce the polarised neutron scattering cross sections,'...
            ' before binning with energy transfer!'])
    end
else
    Pab = [];
end

if parType(5)
    if isfield(spectra,'intP') && ~isempty(spectra.intP)
        intP = spectra.intP;
    else
        error('sw_egrid:WrongInput',['Reference to non-existent field ''intP'','...
            ' use ''sw_neutron(''pol'',true)'' to produce the polarised neutron scattering cross sections,'...
            ' before binning with energy transfer!'])
    end
else
    intP = [];
end

% pack all cross section into a cell for easier looping

if isfield(spectra,'Sab') && iscell(spectra.Sab)
    nTwin = numel(spectra.Sab);
    if ~isfield(spectra,'omega')
        omega = {};
    else
        omega = spectra.omega;
    end
    
    Sab   = spectra.Sab;
%     intP  = {intP};
    if ~iscell(intP); intP = {intP}; end
    if ~iscell(Pab); Pab = {Pab}; end
    if ~iscell(Mab); Mab = {Mab}; end
    if ~iscell(Sperp); Sperp = {Sperp}; end
%     Pab   = {Pab};
%     Mab   = {Mab};
%     Sperp = {Sperp};

elseif isfield(spectra,'Sperp') && iscell(spectra.Sperp)
    nTwin = numel(spectra.Sperp);
    if ~isfield(spectra,'omega')
        omega = {};
    else
        omega = spectra.omega;
    end
    Sab   = repmat({[]},1,nTwin);
    intP  = {intP};
    Pab   = {Pab};
    Mab   = {Mab};
    Sperp = spectra.Sperp;
elseif ~isfield(spectra,'Sab')
    nTwin = 1;
    if ~isfield(spectra,'omega')
        omega = {[]};
    else
        omega = {spectra.omega};
    end
    Sab   = {[]};
    intP  = {intP};
    Pab   = {Pab};
    Mab   = {Mab};
    Sperp = {spectra.Sperp};
else
    nTwin = 1;
    if ~isfield(spectra,'omega')
        omega = {[]};
    else
        omega = {spectra.omega};
    end
    Sab   = {spectra.Sab};
    intP  = {intP};
    Pab   = {Pab};
    Mab   = {Mab};
    Sperp = {Sperp};
end

% number of modes and Q points
if iscell(spectra.Sab)
    nMode = size(spectra.Sab{1},3);
else
    nMode = size(spectra.Sab,3);
end
nHkl  = numel(spectra.hkl)/3;
sHkl  = [size(spectra.hkl) 1];

for ii = 1:nTwin
    Sab{ii}   = reshape(Sab{ii},3,3,nMode,[]);
    if ~isempty(omega{1})
        omega{ii} = reshape(omega{ii},nMode,[]); %#ok<AGROW>
    end
    
    if ~isempty(intP{1})
        intP{ii}   = reshape(intP{ii},3,nMode,[]);
    end
    if ~isempty(Pab{1})
        Pab{ii}   = reshape(Pab{ii},3,3,nMode,[]);
    end
    if ~isempty(Mab{1})
        Mab{ii}   = reshape(Mab{ii},3,3,nMode,[]);
    end
    if ~isempty(Sperp{1})
        Sperp{ii}   = reshape(Sperp{ii},nMode,[]);
    end
end

% Default value of the modeIdx vector selects all modes for output
if isempty(param.modeIdx)
    param.modeIdx = 1:nMode;
end

% DSF stores the intensity that is going to be convoluted
DSF = cell(nConv,nTwin);

for tt = 1:nTwin
    for ii = 1:numel(parsed)
        par0 = parsed{ii};
        DSF{ii,tt} = zeros(nMode, nHkl);
        for jj = 1:length(par0.type)
            switch par0.type{jj}(1)
                case 1
                    DSF{ii,tt} = DSF{ii,tt} + par0.preFact(jj)*Sperp{tt};
                case 2
                    DSF{ii,tt} = DSF{ii,tt} + par0.preFact(jj)*permute(Sab{tt}(par0.type{jj}(2),par0.type{jj}(3),:,:),[3 4 1 2]);
                case 3
                    DSF{ii,tt} = DSF{ii,tt} + par0.preFact(jj)*permute(Mab{tt}(par0.type{jj}(2),par0.type{jj}(3),:,:),[3 4 1 2]);
                case 4
                    DSF{ii,tt} = DSF{ii,tt} + par0.preFact(jj)*permute(Pab{tt}(par0.type{jj}(2),par0.type{jj}(3),:,:),[3 4 1 2]);
                case 5
                    DSF{ii,tt} = DSF{ii,tt} + par0.preFact(jj)*permute(intP{tt}(par0.type{jj}(2),:,:),[2 3 1]);
                otherwise
                    DSF{ii,tt} = DSF{ii,tt} + par0.preFact(jj)*spectra.(par0.type{jj});
            end
        end
    end
end

% save the edge bins
param.Evect = sort(param.Evect);
if eBin
    ebin_edges = param.Evect;
    ebin_cens = (ebin_edges(2:end)+ebin_edges(1:(end-1)))/2;
    dE = diff(ebin_cens);
else
    ebin_cens = param.Evect;
    dE = diff(ebin_cens);
    ebin_edges = [ebin_cens-0.5*[dE(1) dE], ebin_cens(end) + dE(end)/2];
end
nE = numel(ebin_cens);

if isfield(spectra,'omega')   
    if param.imagChk
        % find the maximum of the imaginary part of the spin wave energies
        % checks only the first twin!
        ioMax = max(abs(imag(omega{1}(:))));
        if ioMax > max(abs(dE(:)))
            error('egrid:BadSolution',['The imaginary part of the spin '...
                'wave energes is larger than the bin size! Improve '...
                'your calculation or disable imagChk option!']);
        elseif param.autoEmin
            if abs(ebin_edges(1)) < ioMax
                ebin_edges(1) = ebin_edges(1)+ioMax;
                ebin_cens(1) = ebin_cens(1)+ioMax/2;
            end
        end
    end    
    
    % Calculate Bose temperature factor for magnons
    if param.T==0
        nBose = double(ebin_cens>=0);
    else
        nBose = 1./(exp(abs(ebin_cens)./(spectra.obj.unit.kB*param.T))-1)+double(ebin_cens>=0);
    end
    
    % Sums up the intensities in DSF into swConv.
    swConv = cell(nConv,nTwin);
    
    for tt = 1:nTwin
        for ii = 1:nConv
            real_eigvals = real(omega{tt}(param.modeIdx, :));
            % find energy bin (cen) index coinciding with evals in omega
            ien = discretize(real_eigvals, ebin_edges);
            % set eigvals < zeroEnergyTol to naN (will be ignored)
            izero_eigval = abs(real_eigvals(:)) < param.zeroEnergyTol;
            ien(izero_eigval) = NaN;
             % NaN in ien implies eigvals not in extent of Evect
            ien_valid =  ~isnan(ien(:)); 
            % get hkl index of each ien bin (column index in real_eigvals)
            [~, ihkl] = ind2sub(size(ien), 1:numel(ien));
            % get index of bins in final sxConv field
            % normally ien(ien_valid) is a colulmn vector but not in case
            % of one modeIDx specified so need to reshape
            sw_conv_idx = [reshape(ien(ien_valid), [], 1), ihkl(ien_valid)'];
            if ~isempty(sw_conv_idx)
                % sum intensities and pad energies above max eigval with 0
                DSF_valid = DSF{ii,tt}(param.modeIdx, :);
                DSF_valid(DSF_valid > param.maxDSF) = 0;
                swConv{ii,tt} = accumarray(sw_conv_idx, DSF_valid(ien_valid), [nE, nHkl]);
                % Multiply the intensities with the Bose factor.
                swConv{ii,tt} = bsxfun(@times,swConv{ii,tt},nBose');
                swConv{ii,tt}(isnan(swConv{ii,tt})) = 0;
            else
                swConv{ii,tt} = zeros(numel(ebin_cens), nHkl);
            end
        end
    end
    
else
    swConv = DSF;
    nE     = 1;
end

% sum up twin spectra if requested
if param.sumtwin
    % normalised volume fractions of the twins
    if isfield(spectra,'obj')
        vol = spectra.obj.twin.vol/sum(spectra.obj.twin.vol);
    else
        vol = 1;
    end
    swConvT = cell(nConv,1);
    DSFT    = cell(nConv,1);
    for ii = 1:nConv
        swConvT{ii} = swConv{ii,1}*vol(1);
        DSFT{ii}    = DSF{ii,1}*vol(1);
        for tt = 2:nTwin
            swConvT{ii} = swConvT{ii} + swConv{ii,tt}*vol(tt);
            DSFT{ii}    = DSFT{ii} + DSF{ii,tt}*vol(tt);
        end
    end
    swConv = swConvT;
    DSF    = DSFT;
end


spectra.T     = param.T;
spectra.Evect = ebin_edges;

if isfield(spectra,'swRaw')
    spectra = rmfield(spectra,'swRaw');
end

if numel(swConv) == 1
    spectra.swConv    = reshape(swConv{1},[nE sHkl(2:end)]);
    spectra.swInt     = reshape(DSF{1},[nMode sHkl(2:end)]);
    spectra.component = param.component{1};
else
    spectra.swConv = cell(size(swConv));
    spectra.swInt  = cell(size(DSF));
    for ii = 1:nConv
        spectra.swConv{ii} = reshape(swConv{ii},[nE sHkl(2:end)]);
        spectra.swInt{ii}  = reshape(DSF{ii},[nMode sHkl(2:end)]);
    end
    spectra.component = param.component;
end

spectra.param.sumtwin  = param.sumtwin;

end
