function issym = sw_issymspec(spectra)
% Checks if a given spectrum structure is symbolic or not
%
% ### Syntax
%
% `issym = sw_issymspec(spectra)`
%
% ### Description
%
% Calling spinwave() for a symbolic spinw object results in
% a symbolic spectrum which cannot be used with sw_neutron etc.
% This function checks if a given spectrum is symbolic or not.
%
% ### Input Arguments
%
% `spectra`
% : Input spectra structure.
%

if nargin < 1
    error('sw_issymspec:WrongInput', 'Missing input spectra structure');
end

% Spectrum must have `ham` or `omega` field
if ~isfield(spectra, 'ham')
    if isfield(spectra, 'omega')
        symobj = spectra.omega;
    else
        error('sw_issymspec:WrongInput', 'Invalid input spectra structure');
    end
else
    symobj = spectra.ham;
end
issym = isa(symobj, 'sym');

end
