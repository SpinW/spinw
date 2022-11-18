function spectra = sw_omegasum(spectra, varargin)
% removes degenerate and ghost magnon modes from spectrum
% 
% ### Syntax
% 
% `spec = sw_omegasum(spec,Name,Value)`
% 
% ### Description
% 
% `spec = sw_omegasum(spec,Name,Value)` removes the degenerate modes from
% the dispersion stored in `spec.omega` and sorts the modes according to
% increasing energy. It also removes ghost modes if a lower intensity limit
% is given.
%  
% The degenerate energies are substituted with `NaN` values.
%
% {{warning Be carefull, after this function [sw_egrid] won't work properly.
% This function won't work with spectra of multiple twins.}}
% 
% ### Name-Value Pair Arguments
% 
% `'tol'`
% : Energy tolerance, within the given value two energies are considered
%   equal. Default value is $10^{-5}$.
% 
% `'zeroint'`
% : The minimum intensity value, below which the mode is removed. Default
%   value is 0 (no modes are dropped due to weak intensity).
% 
% `'emptyval'`
% : Value that is assigned to modes that are removed. Default value is NaN
%   (good for plotting). 0 can be used if further numerical analysis, such
%   as binning will be applied.
% 
% ### See Also
% 
% [spinw.spinwave] \| [sw_egrid]
%

if sw_issymspec(spectra)
    error('sw_omegasum:SymbolicInput', 'This function does not handle symbolic spectra');
end

if iscell(spectra.omega)
    error('sw_omegasum:NoTwin','The sw_omegasum() function doesn''t work for spectra calculated for multiple twins!');
end

inpForm.fname  = {'tol' 'zeroint' 'emptyval'};
inpForm.defval = {1e-5  0         NaN       };
inpForm.size   = {[1 1] [1 1]     [1 1]     };

param = sw_readparam(inpForm, varargin{:});

tol = param.tol;
omega    = real(spectra.omega);
omega(isnan(omega)) = 0;
omega(omega<0) = 0;
omegaCol = omega*0+param.emptyval(1);
intCol   = omegaCol;

nQ = size(omega,2);

if isfield(spectra,'swInt') && param.zeroint > 0
    omega(abs(real(spectra.swInt)) < param.zeroint) = 0;
end

for ii = 1:nQ
    oTemp = sw_uniquetol(omega(:,ii)',param.tol);
    oTemp(oTemp == 0) = NaN;
    omegaCol(1:numel(oTemp),ii) = sort(oTemp);
    oTemp(isnan(oTemp)) = [];
    
    if isfield(spectra,'swInt') && ~isempty(oTemp)
        for jj = 1:numel(oTemp)
            intCol(jj,ii) = sum(spectra.swInt(abs(omegaCol(jj,ii)-omega(:,ii))<tol,ii));
        end
    end
    
end

spectra.omega = omegaCol;
spectra.swInt = intCol;
end
