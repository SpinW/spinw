function result = optmagk(obj,varargin)
% determines the magnetic propagation vector
% 
% ### Syntax
% 
% `res = optmagk(obj,Name,Value)`
% 
% ### Description
% 
% `res = optmagk(obj,Name,Value)` determines the optimal propagation vector
% using the Luttinger-Tisza method. It calculates the Fourier transform of
% the Hamiltonian as a function of wave vector and finds the wave vector
% that corresponds to the smalles global eigenvalue of the Hamiltonian.
% The global optimization is achieved using Particle-Swarm optimizer. This
% function sets k and F in spinw.mag_str, and also returns them.
% 
% ### Input Arguments
% 
% `obj`
% : [spinw] object.
% 
% ### Name-Value Pair Arguments
%
% `kbase`
% : Provides a set of vectors that span the space for possible propagation
%   vectors:
%
%   $ \mathbf{k} = \sum_i C(i)\cdot \mathbf{k}_{base}(i);$
%
%   where the optimiser determines the $C(i)$ values that correspond
%      to the lowest ground state energy. $\mathbf{k}_{base}$ is a
%      matrix with dimensions $[3\times n_{base}]$, where $n_{base}\leq 3$. The basis
%      vectors have to be linearly independent.
% 
% The function also accepts all options of [ndbase.pso].
% 
% ### Output Arguments
% 
% `res`
% : Structure with the following fields:
%   * `k`       Value of the optimal k-vector, with values between 0
%                       and 1/2.
%   * `F`       Fourier components for every spin in the magnetic
%                       cell.
%   * `E`       The most negative eigenvalue at the given propagation
%                       vector.
%   * `stat`    Full output of the [ndbase.pso] optimizer.
% 
% ### See Also
% 
% [ndbase.pso]
%

inpForm.fname  = {'kbase'};
inpForm.defval = {[]     };
inpForm.size   = {[3 -1] };
inpForm.soft   = {true   };

warning('off','sw_readparam:UnreadInput')
param = sw_readparam(inpForm, varargin{:});
softparamcheck(["kbase"], 'optmagk', param, varargin{:});
warning('on','sw_readparam:UnreadInput')


% calculate symbolic Fourier transformation if obj is in symbolic mode
if obj.symbolic
    warning('spinw:optmagk:NoSymbolic','The function does not work in symbolic mode!');
    return
end


% generate exchange couplings
[SS, ~, RR] = obj.intmatrix('fitmode',2,'extend',false,'conjugate',true,'zeroC',false);

% add dipolar interactions
SS.all = [SS.all SS.dip];

% list of magnetic atoms in the unit cell
matom = obj.matom;
% number of magnetic atoms in the unit cell
nMagAtom = numel(matom.idx);
% number of bonds
nBond   = size(SS.all,2);

% interacting atom1
%atom1 = SS.all(4,:);
atom1 = mod(SS.all(4,:)-1,numel(matom.S))+1;
% interacting atom1
%atom2 = SS.all(5,:);
atom2 = mod(SS.all(5,:)-1,numel(matom.S))+1;

% exchange energy: J_ij*S_i*S_j
JJ = bsxfun(@times,SS.all(6:14,:),matom.S(atom1).*matom.S(atom2));

% distance vector between interacting atoms in l.u.
dR = SS.all(1:3,:)+RR(:,atom2)-RR(:,atom1);

% number of Q points
nHkl    = 1;

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

if isempty(param.kbase)
    % find the dimensionality of the bonds
    % bond vectors of non-zero couplings
    dl = SS.all(1:3,:);
    % system dimensionality
    D = rank(dl);
    % k-vector directions
    kbase = orth(dl);
else
    kbase = param.kbase;
    D = size(kbase,2);
end

kones = ones(1,D);

warning('off','sw_readparam:UnreadInput')
% optimise the energy using particle swarm
[pOpt0, ~, ~] = ndbase.pso([],@optfun,1/4*kones,'lb',0*kones,'ub',kones,...
    'TolFun',1e-5,'TolX',1e-5,'MaxIter',1e3, varargin{:});
warning('on','sw_readparam:UnreadInput')

% generate an R-value
optfun2 = @(p)(1e-2+optfun(p)-optfun(pOpt0))^2;

% optimize further using LM
[pOpt, ~, stat] = ndbase.lm2([],optfun2,pOpt0,'TolFun',1e-16,'TolX',1e-6);

kOpt = (kbase*pOpt(:));

% is there any value larger than 1/2
kOpt = mod(kOpt,1);

[Eopt, V] = optfun(pOpt);
% 
% % sum up on all atoms
% if isreal(V)
%     n = [];
% else
%     %V = sum(reshape(V,3,[]),2);
%     V = reshape(V,3,[]);
%     % find normal vector
%     %n = cross(real(V(1:3)),imag(V(1:3)));
%     n = cross(real(V),imag(V));
%     n = bsxfun(@rdivide,n,sqrt(sum(n.^2,1)));
%     %n = n(:)';
%     
%     if any(isnan(n))
%         warning('spinw:optmagk:NormalVector','The normal vector is undefined, using [001]!')
%         n = [];
%     end
% end
% 
% if isempty(n)
%     n = repmat([0;0;1],1,nMagAtom);
% end
% 
% % direction of the imaginary moment
% iM = repmat([0;1;0],1,nMagAtom);
% rM = cross(iM,n);
% nZero = sum(~any(rM,1));
% iM(:,~any(rM,1)) = repmat([1;0;0],1,nZero);
% rM = cross(iM,n);

F = reshape(V,3,[]);
% save the optimized values
if isreal(F)
    F = bsxfun(@rdivide,F,sqrt(sum(F.^2,1)));
    % find normal vectors
    v = repmat([0;1;0],1,nMagAtom);
    iM = cross(v,F);
    nZero = sum(~any(iM,1));
    v(:,~any(iM,1)) = repmat([1;0;0],1,nZero);
    iM = cross(v,F);
    iM = bsxfun(@rdivide,iM,sqrt(sum(iM.^2,1)));
    F  = F+1i*iM;
end
obj.genmagstr('mode','fourier','k',kOpt','S',F);

% output results
result.k = kOpt;
result.E = Eopt;
result.F = F;
result.stat = stat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fit function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [E, V] = optfun(p)
        %
        
        kvect = (kbase*p(:))';
        
        % exponents, dimension: 1 x nHkl x nBond
        ExpF = permute(exp(1i*2*pi*sum(bsxfun(@times,kvect,permute(dR,[3 1 2])),2)),[2 1 3]);
        
        % J*exp(ikr), dimensions 9 x nBond x nHkl
        Jexp = permute(bsxfun(@times,permute(JJ,[1 3 2]),ExpF),[1 3 2]);
        
        % Exchange energy: A_ij = S_i*S_j*J_ij*exp(ik*2pi*(R_i-R_j))
        % Summed up on all bonds
        % Dimensions: 9 x nMagExt x nMagExt x nHkl
        %                 atom1     atom2
        ft = accumarray(idxAll,Jexp(:)',[9 nMagAtom nMagAtom nHkl]);
        
        if nargout == 1
            % reshape into 3 x 3 x nMagExt x nMagExt x nHkl
            E = min(real(eig(reshape(permute(reshape(ft,[3 3 nMagAtom nMagAtom nHkl]),[1 3 2 4 5]),3*nMagAtom,3*nMagAtom,[]))));
            
        else
            % reshape into 3 x 3 x nMagExt x nMagExt x nHkl
            [V,E] = eig(reshape(permute(reshape(ft,[3 3 nMagAtom nMagAtom nHkl]),[1 3 2 4 5]),3*nMagAtom,3*nMagAtom,[]));
            % find degenerate eigenvalues
            E = diag(E);
            [E,idxE] = sort(real(E));
            V = V(:,idxE);
            sel = (E-E(1))<10*eps;
            V = sum(V(:,sel),2);
            E = E(1);
        end
        
        E = real(E);
        
    end

end