function res = fourier(obj,hkl,varargin)
% calculates the Fourier transformation of the Hamiltonian
% 
% ### Syntax
% 
% `F = fourier(obj,Q,Name,Value)`
% 
% ### Description
% 
% `F = fourier(obj,hkl,Name,Value)` calculates the following Fourier sum:
%
% $J(\mathbf{k}) = \sum_{i,j} J_{i,j} * \exp(i \mathbf{k}\cdot \mathbf{d}_{i,j})$
%
% The code is optimised for calculating the sum for large number of wave
% vectors and alternatively for a large number of $d_{i,j}$ vectors (large
% system size). The single ion anisotropy is not included in the sum.
% 
% ### Input Arguments
% 
% `obj`
% : [spinw] object.
% 
% `Q`
% : Defines the $Q$ points where the spectra is calculated, in reciprocal
%   lattice units, size is $[3\times n_{Q}]$. $Q$ can be also defined by
%   several linear scan in reciprocal space. In this case `Q` is cell type,
%   where each element of the cell defines a point in $Q$ space. Linear scans
%   are assumed between consecutive points. Also the number of $Q$ points can
%   be specified as a last element, it is 100 by defaults. 
%   
%   For example to define a scan along $(h,0,0)$ from $h=0$ to $h=1$ using
%   200 $Q$ points the following input should be used:
%   ```
%   Q = {[0 0 0] [1 0 0]  50}
%   ```
%
%   For symbolic calculation at a general reciprocal space point use `sym`
%   type input. 
%
%   For example to calculate the spectrum along $(h,0,0)$ use:
%   ```
%   Q = [sym('h') 0 0]
%   ```
%   To calculate spectrum at a specific $Q$ point symbolically, e.g. at
%   $(0,1,0)$ use:
%   ```
%   Q = sym([0 1 0])
%   ```
% 
% ### Name-Value Pair Arguments
% 
% `'extend'`
% : If `true`, the Fourier transform will be calculated on the
%   magnetic supercell, if `false` the crystallographic cell will
%   be considered. Default is `true.`
% 
% `'isomode'`
% : Defines how Heisenberg/non-Heisenberg Hamiltonians are
%   treated. Can have the following values:
%   * `'off'`   Always output the $[3\times 3]$ form of the
%               Hamiltonian, (default).
%   * `'auto'`  If the Hamiltonian is Heisenberg, only output
%               one of the diagonal values from the $[3\times 3]$
%               matrices to reduce memory consumption.
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
% `res` struct type with the following fields:
% * `ft`        contains the Fourier transform in a matrix with dimensions
%               $[3\times 3\times n_{magExt}\times n_{magExt}\times
%               n_{hkl}]$ or $[1\times 1\times n_{magExt}\times n_{magExt}\times n_{hkl}]$
%               for Heisenberg and non-Heisenberg Hamiltonians respectively
%               (if isomode is `'auto'`). Here $n_{magExt}$ is the number of
%               magnetic atoms in the magnetic cell and $n_{hkl}$ is the number
%               of reciprocal space points.
% * `hkl`       Matrix with the given reciprocal space points stored in a
%               matrix with dimensions $[3\times n_{hkl}]$.
% * `isiso`     True is the output is in Heisenberg mode, when the `ft`
%               matrix has dimensions of $[1\times 1\times n_{magExt}\times n_{magExt}\times n_{hkl}]$,
%               otherwise it is `false`.
% 
% ### See Also
% 
% [spinw.optmagk]
%

% TODO: test for magnetic supercell

inpForm.fname  = {'extend' 'fid' 'isomode' 'sublat'};
inpForm.defval = {true     -1    'off'     []      };
inpForm.size   = {[1 1]    [1 1] [1 -1]    [1 -2]  };
inpForm.soft   = {false    false false     true    };

param = sw_readparam(inpForm, varargin{:});
pref = swpref;

switch param.isomode
    case 'auto'
        isomode = true;
    case 'off'
        isomode = false;
    otherwise
        error('spinw:WrongInput','The given isomode option is invalid!')
end

if param.fid == -1
    % Print output into the following file
    fid = pref.fid;
else
    fid = param.fid;
end

% for linear scans create the Q line(s)
if nargin > 1
    if iscell(hkl)
        hkl = sw_qscan(hkl);
    elseif numel(hkl)==3
        hkl = hkl(:);
    end
else
    hkl = [];
end

% calculate symbolic Fourier transformation if obj is in symbolic mode
if obj.symbolic
    if numel(hkl) == 3
        hkl = sym(hkl);
    end
    
    if ~isa(hkl,'sym')
        inpForm.fname  = {'fitmode'};
        inpForm.defval = {false    };
        inpForm.size   = {[1 1]    };
        param0 = sw_readparam(inpForm, varargin{:});
        
        if ~param0.fitmode
            fprintf0(fid,['No symbolic hkl value was given, Fourier'...
                ' transformation for general Q (h,k,l) will be calculated!\n']);
        end
        res = obj.fouriersym(varargin{:});
    else
        res = obj.fouriersym(varargin{:},'hkl',hkl);
    end
    return
end

% generate exchange couplings
%[SS, ~, RR] = obj.intmatrix('fitmode',true,'extend',param.extend,'conjugate',true);
[SS, ~, RR] = obj.intmatrix(struct('fitmode',true,'extend',param.extend,'conjugate',true),'noCheck');

% is there Heisenberg exchange only
matidx = unique(obj.coupling.mat_idx(:)');
matidx(matidx==0) = [];
isIso = all(sw_mattype(obj.matrix.mat(:,:,matidx))==1);

fprintf0(fid,'Calculating the Fourier transformation of the spin Hamiltonian...\n')
if isIso
    fprintf0(fid,'Heisenberg model...\n')
else
    fprintf0(fid,'Non-Heisenberg model...\n')
end
% list of magnetic atoms in the unit cell
matom = obj.matom;
% number of Q points
nHkl    = size(hkl,2);
if param.extend
    nExt = double(obj.mag_str.nExt);
else
    nExt = 1;
end
% number of magnetic atoms in the crystallographic unit cell
nMag    = numel(matom.idx);
% number of bonds
nBond   = size(SS.all,2);

% interacting atom1
atom1 = mod(SS.all(4,:)-1,nMag)+1;
% interacting atom2
atom2 = mod(SS.all(5,:)-1,nMag)+1;

% magnetic couplings, 3x3xnJ
%JJ = cat(3,reshape(SS.all(6:14,:),3,3,[]),SI.aniso);

if isIso
    % keep only a single value for Heisenberg interactions
    JJ = SS.all(6,:);
else
    JJ = SS.all(6:14,:);
end

% include the spin values in the exchange energy: J_ij*S_i*S_j
JJ = bsxfun(@times,JJ,matom.S(atom1).*matom.S(atom2));

% distance vector between interacting atoms in l.u.
dR = SS.all(1:3,:)+RR(:,atom2)-RR(:,atom1);

% exponents, dimension: 1 x nHkl x nBond
ExpF = permute(exp(1i*2*pi*sum(bsxfun(@times,hkl',permute(dR,[3 1 2])),2)),[2 1 3]);

% J*exp(ikr), dimensions 9 x nBond x nHkl (1 x nBond x nHkl for Heisenberg)
Jexp = permute(bsxfun(@times,permute(JJ,[1 3 2]),ExpF),[1 3 2]);

% reduce the number of sublattices
if ~isempty(param.sublat)
    atom1 = param.sublat(atom1);
    atom2 = param.sublat(atom2);
    nMag  = max(param.sublat);
    nMag0 = numel(obj.matom.idx);
    nsubl = nMag / nMag0;
    fprintf0(fid,'Remapping magnetic atoms into a new set of sublattices...\n');
else
    nsubl = 1;
end

% number of magnetic atoms in the magnetic supercell
nMagExt = prod(nExt)*nMag;

% optimize for large number of hkl, is not faster for more than ~10
% bonds% if nHkl>nBond
%     ft = zeros(3,3,nMagExt,nMagExt,nHkl);
%     if isIso
%         for ii = 1:nBond
%             ft(1,1,atom1(ii),atom2(ii),:) = ft(1,1,atom1(ii),atom2(ii),:)+permute(Jexp(1,ii,:),[1 2 4 5 3]);
%         end
%         ft(2,2,:,:,:) = ft(1,1,:,:,:);
%         ft(3,3,:,:,:) = ft(1,1,:,:,:);
%     else
%         for ii = 1:nBond
%             ft(:,:,atom1(ii),atom2(ii),:) = reshape(Jexp(:,ii,:),3,3,[]);
%         end
%     end
% end

% optimize for large number of bonds
if ~isIso
    % general, non-Heisenberg model
    % J^\alpha^\beta component indices in numbers 1:9
    idx1 = repmat((1:9)',[1 nBond nHkl]);
    % atom1 indices
    idx2 = repmat(atom1,[9 1 nHkl]);
    % atom2 indices
    idx3 = repmat(atom2,[9 1 nHkl]);
    % hkl indices
    idx4 = repmat(permute(1:nHkl,[1 3 2]),[9 nBond 1]);
    % indices all
    idxAll = [idx1(:) idx2(:) idx3(:) idx4(:)];
    
    % Exchange energy: A_ij = S_i*S_j*J_ij*exp(ik*2pi*(R_i-R_j))
    % Summed up on all bonds
    % Dimensions: 9 x nMagExt x nMagExt x nHkl (1 x nMagExt x nMagExt x nHkl for Heisenberg)
    %                 atom1     atom2
    % the result is always real, because we use the conjugate of the bonds as
    % well
    ft = real(accumarray(idxAll,Jexp(:)',[9 nMagExt nMagExt nHkl]));
    % reshape into 3 x 3 x nMagExt x nMagExt x nHkl
    ft = reshape(ft,[3 3 nMagExt nMagExt nHkl]);
else
    % Heisenberg model
    % J^\alpha^\beta component indices is always 1
    % ...
    % atom1 indices
    idx2 = repmat(atom1,[1 1 nHkl]);
    % atom2 indices
    idx3 = repmat(atom2,[1 1 nHkl]);
    % hkl indices
    idx4 = repmat(permute(1:nHkl,[1 3 2]),[1 nBond 1]);
    % indices all
    idxAll = [idx2(:) idx3(:) idx4(:)];
    
    % Exchange energy: A_ij = S_i*S_j*J_ij*exp(ik*2pi*(R_i-R_j))
    % Summed up on all bonds
    % Dimensions: 9 x nMagExt x nMagExt x nHkl (1 x nMagExt x nMagExt x nHkl for Heisenberg)
    %                 atom1     atom2
    % the result is always real, because we use the conjugate of the bonds as
    % well
    ft = real(accumarray(idxAll,Jexp(:)',[nMagExt nMagExt nHkl]));
    % reshape into 3 x 3 x nMagExt x nMagExt x nHkl
    ft = reshape(ft,[1 1 nMagExt nMagExt nHkl]);
    if ~isomode
        ft(3,3,:,:,:) = ft(1,1,:,:,:);
        ft(2,2,:,:,:) = ft(1,1,:,:,:);
    end
end

% save results in a struct
% scale ft with the number of sublattices
res.ft    = ft * nsubl;
res.hkl   = hkl;
% Heisenberg output
res.isiso = size(ft,1)==1;

fprintf0(fid,'Calculation finished.\n');

end
