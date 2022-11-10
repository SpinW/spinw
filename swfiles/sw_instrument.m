function spectra = sw_instrument(spectra, varargin)
% convolutes spectrum with resolution function
% 
% ### Syntax
% 
% `spectra = sw_instrument(spectra,Name,Value)`
% 
% ### Description
% 
% `spectra = sw_instrument(spectra,Name,Value)` can convolute an energy
% binned spectrum with different energy resolution functions and add other
% effects that introduced by measurement (such as the kinematic limit for
% neutron scattering, finite momentum resolution or finite detector
% coverage).
%  
% 
% ### Name-Value Pair Arguments
% 
% `'dE'`
% : Convolutes the spectrum with a Gaussian in energy, where the width is
%   defined by the FWHM value. The accepted values are:
%   * *string*   File name, that contains the FWHM energy
%                resolution values as a function of energy
%                transfer. The file has to contain two columns,
%                first is the energy values, the second is the
%                FWHM resolution at the given energy transfer
%                value, see [sw_res] function for details.
%   * *number*   Constant FWHM energy resolution given by the number.
%   * *matrix*   Dimensions of $[N\times 2]$, where the first column contains the
%                energy transfer values, second column contains
%                the FWHM resolution values. These discrete values will
%                be fitted using a polynomial with a fixed
%                degree, see [sw_res] for details.
%   * *function* Function handle of a resolution function
%   with the following header:
%   ```
%   E_FWHM = res_fun(E)
%   ```
%   where `E_FWHM` is the FWHM energy resolution and `E` is the energy transfer
%   value.
% 
% `'func'`
% : Shape of the energy resolution function if different from Gaussian.
%   For details see [sw_resconv].
% 
% `'polDeg'`
% : Degree of the polynomial that is fitted to the discrete energy 
%   resolution data. Only used if `dE` is a matrix of string. Default value
%   is 5.
% 
% `'dQ'`
% : Momentum transfer resolution of the instrument, FWHM is
%   given in \\ang$^{-1}$ units by default, unless different units
%   are defined in [spinw.unit]. Default value is 0 for no convolution.
% 
% `'thetaMin'`
% : Minimum scattering angle in \\deg, default value is 0. Can be only
%   applied if one of the `ki`, `Ei`, `kf` or `Ef` parameters is defined.
% 
% `'plot'`
% : If the resolution is read from file and plot option is
%   true, the energy dependent resolution values together with the
%   polynomial fit will be plotted in a new figure. Default value is
%   `true`.
% 
% `'norm'`
% : If true, the data is normalized to mbarn units. Default is
%   false. If no g-tensor is included in the spin wave
%   calculation, $g = 2$ will be assumed for the conversion.
% 
% `'useRaw'`
% : If `false`, the already modified `spectra.swConv` field is
%   modified further instead of the original powder spectrum
%   stored in `spectra.swRaw`. Default value is `true`.
% 
% For simulating the effect of the neutron kinematic limit or the finite 
% detector coverage of a neutron spectrometer one of the following
% parameter has to be given. The unit of these quantities is defined in
% [spinw.unit] with default momentum unit of \\ang$^{-1}$ and energy
% unit of meV.
%
% `'ki'`
% : Fixed momentum of the incident neutrons.
% 
% `'Ei'`
% : Fixed energy of the incident neutrons.
% 
% `'kf'`
% : Fixed final momentum of the neutrons.
% 
% `'Ef'`
% : Fixed final energy of the neutrons.
% 
% `'fid'`
% : Defines whether to provide text output. The default value is determined
%   by the `fid` preference stored in [swpref]. The possible values are:
%   * `0`   No text output is generated.
%   * `1`   Text output in the MATLAB Command Window.
%   * `fid` File ID provided by the `fopen` command, the output is written
%           into the opened file stream.
%
% ### Output Arguments
% 
% `spectra`
% : Struct variable, same as input with following additional fields:
% * `norm`      `true`, if the spectrum is normalised to mbarn units.
% * `ki`        Fixed incident neutron wave vector if defined in the input.
% * `kf`        Fixed final neutron wave vector if defined in the input.
% * `dE`        Energy resolution polynomial as given in the input.
% * `dQ`        FWHM of the momentum resolution.
% * `swRaw`     Original simulated data before the application of
%               `sw_instrument`.
% 
% ### See Also
% 
% [polyfit] \| [polyval] \| [sw_res] \| [sw_resconv]
%
% *[FWHM]: Full Width at Half Maximum
%

if nargin == 0
    swhelp sw_instrument
    return
end

if sw_issymspec(spectra)
    error('sw_instrument:SymbolicInput', 'This function does not handle symbolic spectra');
end

func0 = @swfunc.gaussfwhm;

inpForm.fname  = {'dE'    'ki'  'Ei'  'kf'  'Ef'  'plot' 'polDeg' 'thetaMin'};
inpForm.defval = {[]      0     0     0     0     false   5        0        };
inpForm.size   = {[-1 -2] [1 1] [1 1] [1 1] [1 1] [1 1]  [1 1]    [1 1]     };
inpForm.soft   = {true    false false false false false  false    false     };

inpForm.fname  = [inpForm.fname  {'formFact' 'dQ'  'norm' 'useRaw' 'func' 'fid'}];
inpForm.defval = [inpForm.defval { 'auto'    0     false   true    func0  -1   }];
inpForm.size   = [inpForm.size   { [1 -2]    [1 1] [1 1]   [1 1]   [1 1]  [1 1]}];
inpForm.soft   = [inpForm.soft   {false      false false   false   false  false}];

param = sw_readparam(inpForm, varargin{:});
pref = swpref;

% Print output
if param.fid == -1
    fid = pref.fid;
else
    fid = param.fid;
end

if isfield(spectra,'swRaw')
    % take raw convoluted spectra if exists and request
    if param.useRaw
        spectra.swConv = spectra.swRaw;
    end
else
    % otherwise saw the spectra into swRaw field for future use
    spectra.swRaw = spectra.swConv;
end

% remove non-finite values
if ~iscell(spectra.swConv)
    spectra.swConv = {spectra.swConv};
end

% number of plots
nPlot = numel(spectra.swConv);

for ii = 1:nPlot
    spectra.swConv{ii}(~isfinite(spectra.swConv{ii})) = 0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% energy resolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% determine whether to include resolution
calcres = true;
try %#ok<TRYNC>
    if isempty(param.dE)
        calcres = false;
    end
end

if isa(param.dE,'function_handle')
    % do nothing
elseif numel(param.dE)>1
    
    % determine the energy resolution curve from a file or the given matrix
    param.dE = sw_res(param.dE,param.polDeg,param.plot);
    
end

% center bin positions
cEvect = (spectra.Evect(1:(end-1))+spectra.Evect(2:end))/2;

% include resolution
if calcres
    for jj = 1:nPlot
        %spectra.swConv{jj} = sw_resconv(spectra.swConv{jj},cEvect',param.dE);
        spectra.swConv{jj} = sw_resconv(spectra.swConv{jj},cEvect',param.dE,param.func);
    end
    
    fprintf0(fid,'Finite instrumental energy resolution is applied.\n');
end
spectra.dE = param.dE;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Q resolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Q = spectra.hklA;

if size(Q,1) > 1
    %Q = sqrt(sum(bsxfun(@minus,Q,Q(:,1)).^2,1));
    Q = sqrt(sum(Q.^2,1));
end

if param.dQ > 0
    Qconv = [Q(1) Q(1)+cumsum(abs(diff(Q)))];
    
    % standard deviation of the Q resolution Gaussian
    stdG = param.dQ/2.35482;
    
    for jj = 1:nPlot
        swConv = spectra.swConv{jj};
        swConvTemp = swConv * 0;
        for ii = 1:numel(Qconv)
            % Gaussian with intensity normalised to 1, centered on E(ii)
            fG = exp(-((Qconv-Qconv(ii))/stdG).^2/2);
            fG = fG/sum(fG);
            swConvTemp = swConvTemp + swConv(:,ii) * fG;
            
        end
        spectra.swConv{jj} = swConvTemp;
    end
    
    fprintf0(fid,'Finite instrumental momentum resolution of %5.3f A-1 is applied.\n',param.dQ);
end

spectra.dQ = param.dQ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% energy transfer range
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (param.ki>0) || (param.Ei>0)
    FX = 1;
    param.k = param.ki;
    param.E = param.Ei;
elseif (param.kf>0) || (param.Ef>0)
    FX = 2;
    param.k = param.kf;
    param.E = param.Ef;
else
    FX = 0;
end

if FX>0 && param.k == 0
    param.k = sw_converter(param.E,'meV','k');
end

% save ki/kf value into the output
switch FX
    case 1
        spectra.ki = param.k;
        kstr = 'ki';
    case 2
        spectra.kf = param.k;
        kstr = 'kf';
end

if FX > 0
    k0 = param.k;
    cosT = cosd(param.thetaMin);
    sinT = sind(param.thetaMin);
    
    
    for jj = 1:nPlot
        switch FX
            case 1
                % fix ki
                Emax = (k0^2-(k0*cosT-sqrt(Q.^2-k0^2*sinT^2)).^2) * sw_converter(1,'k','meV');
                Emin = (k0^2-(k0*cosT+sqrt(Q.^2-k0^2*sinT^2)).^2) * sw_converter(1,'k','meV');
            case 2
                % fix kf
                Emax = -(k0^2-(k0*cosT+sqrt(Q.^2-k0^2*sinT^2)).^2) * sw_converter(1,'k','meV');
                Emin = -(k0^2-(k0*cosT-sqrt(Q.^2-k0^2*sinT^2)).^2) * sw_converter(1,'k','meV');
        end
        %Emax = (ki^2-(ki*cosT-sqrt(Q.^2-ki^2*sinT^2)).^2) * sw_converter(1,'k','meV');
        %Emin = (ki^2-(ki*cosT+sqrt(Q.^2-ki^2*sinT^2)).^2) * sw_converter(1,'k','meV');
        
        Emax(abs(imag(Emax))>0) = 0;
        Emin(abs(imag(Emin))>0) = 0;
        
        Elist = repmat(cEvect',[1 size(spectra.swConv{jj},2)]);
        Emin  = repmat(Emin,[size(spectra.swConv{jj},1) 1]);
        Emax  = repmat(Emax,[size(spectra.swConv{jj},1) 1]);
        
        
        idx = (Elist<Emin) | (Elist>Emax);
        swConv = spectra.swConv{jj};
        swConv(idx) = NaN;
        spectra.swConv{jj} = swConv;
    end
    fprintf0(fid,'Energy transfer is limited to instrument, using %s=%5.3f A-1.\n',kstr,k0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normalise spectrum to mbarn/meV/f.u. or mbarn/meV/cell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(spectra,'gtensor')
    spectra.gtensor = false;
end

if param.norm
    % Lande' g-factor
    if spectra.gtensor
        % g-tensor is already included in the spinwave calculation
        g = 1;
    else
        % use simple g=2 here
        g = 2;
    end
    % constant: p = gamma*r0/2
    % neutron magnetic moment constant: M = gamma*gammaN
    gamma = 1.91304272; % 1/s/T
    % classical radius of the electron
    r0 = 2.8179403267e-15; % m
    % cross section constant in mbarn
    p2 = (g*gamma*r0/2)^2*1e28*1e3; % mbarn
    
    % convert intensity to mbarn/meV units using the energy bin size
    dE = diff(spectra.Evect);
    % the new Evect is in edge bin mode
    %dE = [dE(1) (dE(2:end)+dE(1:end-1))/2 dE(end)];
    
    for jj = 1:nPlot
        spectra.swConv{jj} = spectra.swConv{jj}*p2./repmat(dE',[1 size(spectra.swConv{jj},2)]);
    end
    
    % set 'normalized units' switch on
    spectra.norm = true;
    if spectra.obj.unit.nformula > 0
        fprintf0(fid,'Intensity is converted to mbarn/meV/f.u. units.\n');
    else
        fprintf0(fid,'Intensity is converted to mbarn/meV/cell units.\n');
    end
    if spectra.gtensor
        fprintf0(fid,'g-tensor was already included in the spin wave calculation.\n');
    else
        fprintf0(fid,'Isotropic g-tensor of 2 assumed here.\n');
    end
    
else
    spectra.norm = false;
end

if nPlot == 1
    spectra.swConv = spectra.swConv{1};
end

end
