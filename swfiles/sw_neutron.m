function spectra = sw_neutron(spectra, varargin)
% calculates neutron scattering cross section
% 
% ### Syntax
% 
% `spectra = sw_neutron(spectra,Name,Value)`
% 
% ### Description
% 
% `spectra = sw_neutron(spectra,Name,Value)` calculates the neutron
% scattering cross section for polarised and unpolarised neutrons. The
% function reads the calculated spin-spin correlation function
% $\mathcal{S}^{\alpha\beta}(\mathbf{Q},\omega)$ and calculates the neutron
% scattering cross section for unpolarized neutrons using the formula:
%
% $S_\perp(Q,\omega)=\sum_{\alpha\beta}(1-\hat{q}^\alpha\hat{q}^\beta)\cdot S^{\alpha\beta}(Q,\omega)$
%  
% It also calculates spin-spin correlation function in the Blume-Maleev
% coordinate system and the complete polarised neutron scattering cross
% section.
% 
% {{note The Blume-Maleev coordinate system is a cartesian coordinate
% system with $x_{BM}$, $y_{BM}$ and $z_{BM}$ basis vectors defined as:
% <br> $x_{BM}$    parallel to the momentum transfer $Q$,
% <br> $y_{BM}$    perpendicular to $x_{BM}$ in the scattering plane,
% <br> $z_{BM}$    perpendicular to the scattering plane.
% }}
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
% `'n'`
% : Normal vector to the scattering plane, in real space ($xyz$
%   coordinate system), stored in a row vector with 3 elements. Default
%   value is `[0 0 1]`.
%
% `'uv'`
% : Cell, that contains two vectors defining the scattering 
%   plane in rlu. If given overwrites the `n` parameter value. For example:
%   `{[1 0 0] [0 1 0]}` stands for the $(h,k,0)$ scattering plane.
% 
% `'pol'`
% : If `true` the cross sections in the Blume-Maleev
%   coordinate system will be also calculated (`inP`, `Pab` and `Mab`
%   fields of the output `spectra`). Default value is `false`.
% 
% ### Output Arguments
% 
% `spectra`
% : Same as the input `spectra` plus the following additional fields:
%   * `param`   Input parameters.
%   * `Sperp`   $S_\perp(i_{mode},\mathbf{Q})$ unpolarised neutron 
%               scattering cross section, stored in a matrix with
%               dimensions of $[n_{mode}\times n_{hkl}]$.
%   * `intP`    $I_P(P_i,i_{mode},\mathbf{Q})$ polarised neutron scattering 
%               cross section, when only the incident neutron polarization
%               is analyzed. It is stored in a matrix with dimensions of
%               $[3\times n_{mode}\times n_{hkl}]$.
%   * `Pab`     $I_{Pab}(P_f,P_i,i_{mode},\mathbf{Q})$ complete polarised 
%               neutron scattering cross section, when the polarisation of
%               both the incident ($P_i$) and the scattered ($P_f$)
%               neutrons are analyzed. Stored in a matrix with dimensions
%               of $[3\times 3\times n_{mode}\times n_{hkl}]$.
%   * `Mab`     $M_{ab}(P_f,P_i,i_{mode},\mathbf{Q})$ components of the 
%               spin-spin correlation function in the blume-Maleev
%               coordinate system, stored in a matrix with dimensions of
%               $[3\times \times3 n_{mode}\times n_{hkl}]$.
%
% If several domains exist in the sample, `Sperp`, `intP`, `Pab` and `Mab`
% will be packaged into a cell, that contains $n_{twin}$ number of
% matrices.
%
% The meaning of the indices above:
% * $P_i$: index of incident polarisation ($1=x_{BM}$, $2=y_{BM}$ or $3=z_{BM}$),
% * $P_f$: index of final polarisation ($1=x_{BM}$, $2=y_{BM}$ or $3=z_{BM}$),
% * $i_{mode}$: index of spin wave mode,
% * $\mathbf{Q}$: index of momentum transfer.
%
% 
% ### See Also
% 
% [sw_egrid] \| [spinw] \| [spinw.spinwave]
%
% *[rlu]: reciprocal lattice units
%

if nargin == 0
    swhelp sw_neutron
    return
end

if sw_issymspec(spectra)
    error('sw_neutron:SymbolicInput', 'This function does not handle symbolic spectra');
end

inpForm.fname  = {'n'     'pol' 'uv'  };
inpForm.defval = {[0 0 1] false  {}    };
inpForm.size   = {[1 3]   [1 1] [1 2] };
inpForm.soft   = {false   false true  };

param = sw_readparam(inpForm,varargin{:});

% normal to the scattering plane in xyz coordinate system
if numel(param.uv) == 0
    n = param.n;
else
    bv = spectra.obj.basisvector;
    u = param.uv{1};
    v = param.uv{2};
    if (numel(u)~=3) || (numel(v)~=3)
        error('sw_neutron:WrongInput','The u and v vectors of the scattering plane are wrong!');
    end
    n = cross(u*2*pi*inv(bv),v*2*pi*inv(bv)); %#ok<MINV>
    n = round(n*1e12)/1e12;
    n = n/norm(n);
end


% input parameters extracted from spectra
% swSpec.Sab   3 x 3 x nMode x nHkl
% swSpec.hklA  3 x nHkl
hklA    = reshape(spectra.hklA,3,[]);
nHkl    = size(hklA,2);
sHkl    = size(spectra.hklA);

% loop over the twins
if ~iscell(spectra.Sab)
    spectra.Sab = {spectra.Sab};
end

nTwin = numel(spectra.Sab);
spectra.intP  = cell(1,nTwin);
spectra.Pab   = cell(1,nTwin);
spectra.Mab   = cell(1,nTwin);
spectra.Sperp = cell(1,nTwin);

for ii = 1:nTwin
    nMode = size(spectra.Sab{ii},3);
    Sab   = reshape(spectra.Sab{ii},[3,3,nMode,nHkl]);
    
    
    % divide Sab into symmetric and anti-symmetric components
    SabA = (Sab - permute(Sab,[2 1 3 4]))/2;
    SabS = (Sab + permute(Sab,[2 1 3 4]))/2;
    
    
    % unplarised nuetron scattering
    
    % Normalized scattering wavevector in xyz coordinate system.
    hklAN = bsxfun(@rdivide,hklA,sqrt(sum(hklA.^2,1)));
    
    % avoid NaN for Q=0
    NaNidx = find(any(isnan(hklAN)));
    for jj = 1:numel(NaNidx)
        if NaNidx(jj) < size(hklAN,2)
            hklAN(:,NaNidx(jj)) = hklAN(:,NaNidx(jj)+1);
        else
            hklAN(:,NaNidx(jj)) = [1;0;0];
        end
    end
    
    hkla = repmat(permute(hklAN,[1 3 2]),[1 3 1]);
    hklb = repmat(permute(hklAN,[3 1 2]),[3 1 1]);
    
    % Perpendicular part of the scattering wavevector.
    qPerp = repmat(eye(3),[1 1 nHkl])- hkla.*hklb;
    qPerp = repmat(permute(qPerp,[1 2 4 3]),[1 1 nMode 1]);
    
    % Dynamical structure factor for neutron scattering
    % Sperp: nMode x nHkl.
    Sperp = permute(sumn(qPerp.*SabS,[1 2]),[3 4 1 2]);
    
    if param.pol
        % polarised neutron scattering (without final polarisation analysis)
        % Levi-Civita symbol
        lcMat = zeros(3,3,3); lcMat([8 12 22]) = 1; lcMat([6 16 20]) = -1;
        
        % polarised scattering in the (xBM,yBM,zBM) coordinate system
        xBM = hklAN;
        zBM = repmat(n',[1 nHkl]);
        
        if norm(zBM(:,1)'*xBM)/nHkl > 1e-5
            warning('sw_neutron:ScatteringPlaneProblem','The normal to the scattering plane is not perpendicular to Q!');
        end
        
        yBM = cross(zBM,xBM);
        yBM = bsxfun(@rdivide,yBM,sqrt(sum(yBM.^2,1)));
        
        % matrix of incoming polarisation
        % Pi: 3 (x,y,z) x nHkl x 3 (xBM,yBM,zBM)
        Pi = cat(3,xBM,yBM,zBM);
        
        % see written notes for the matrix notation!
        % Qdelta: 3(alpha) x 3(gamma) x 3(delta) x nHkl
        Qdelta = repmat(permute(xBM,[3 4 1 2]),[3 3 1 1]);
        
        % mat1: 3(gamma) x 3(alpha) x 3(beta) x nHkl
        mat1 = repmat(permute(sum(repmat(lcMat,[1 1 1 nHkl]).*Qdelta,3),[2 1 3 4]),[1 1 3 1]).*Qdelta;
        % mat2: 3(gamma) x 3(alpha) x 3(beta) x nHkl
        mat2 = repmat(lcMat,[1 1 1 nHkl]);
        
        % mat3: 3 (gamma) x nMode x nHkl
        mat3 = permute(sumn(repmat(permute(2*mat1 + mat2,[2 3 5 4 1]),[1 1 nMode 1 1]).*repmat(SabA,[1 1 1 1 3]),[1 2]),[5 3 4 1 2]);
        Pgamma = repmat(permute(Pi,[1 4 2 3]),[1 nMode 1 1]);
        
        % B:  P * <Sperp x Sperp(t)>
        % B: 3(Pix,Piy,Piz) x nMode x nHkl
        B = permute(sum(Pgamma.*repmat(mat3,[1 1 1 3]),1),[4 2 3 1]);
        % intP: 3(Pix,Piy,Piz) x nMode x nHkl
        intP = 1i*B + repmat(permute(Sperp,[3 1 2]),[3 1 1]);
        
        % final polarisation of the neutrons
        % P*Q: 3(x,y,z) x nHkl x 3 (Pix,Piy,Piz)
        Qbeta = repmat(xBM,[1 1 3]);
        PQ    = repmat(sum(Pi.*Qbeta,1),[3 1 1]);
        % mat4: 3(alpha) x 3(beta) x nHkl x 3 (Pix,Piy,Piz)
        mat4  = repmat(permute(Pi - PQ.*Qbeta,[4 1 2 3]),[3 1 1 1]);
        % mat5: 3(gamma) x 3(beta) x nHkl x 3 (Pix,Piy,Piz)
        mat5  = -mat4.* repmat(permute(Qbeta,[1 4 2 3]),[1 3 1 1]);
        
        % Pab1: 3(x,y,z) x 3(Pix,Piy,Piz) x nMode x nHkl
        Pab1 = permute(sum(repmat(permute(mat4,[1 2 5 3 4]),[1 1 nMode 1 1]).*repmat(SabS,[1 1 1 1 3]),2),[1 5 3 4 2]);
        % Pab2: 3(x,y,z) x 3(Pix,Piy,Piz) x nMode x nHkl
        Pab2 = repmat(permute(sumn(repmat(permute(mat5,[1 2 5 3 4]),[1 1 nMode 1 1]).*repmat(SabS,[1 1 1 1 3]),[1 2]),[1 5 3 4 2]),[3 1 1 1 1]);
        Pab2 = Pab2.*repmat(permute(Qbeta,[1 3 4 2]),[1 1 nMode 1]);
        % Pi*Sperp
        % Pab3: 3(x,y,z) x 3(Pix,Piy,Piz) x nMode x nHkl
        Pab3 = repmat(permute(Pi,[1 3 4 2]),[1 1 nMode 1]).*repmat(permute(Sperp,[3 4 1 2]),[3 3 1 1]);
        % i*<Sperp x Sperp(t)>
        % Pab4: 3(x,y,z) x 3(Pix,Piy,Piz) x nMode x nHkl
        Pab4 = repmat(permute(mat3,[1 4 2 3]),[1 3 1 1]);
        % Pab : 3(x,y,z) x 3(Pix,Piy,Piz) x nMode x nHkl
        Pab = 2*(Pab1 + Pab2) - Pab3 - 1i*Pab4;
        
        % Pab first dimension is in real space --> change to (Pfx,Pfy,Pfz) coordinate
        Pf = permute(Pi,[1 3 2]);
        % invPf converts (x,y,z) --> (Pfx,Pfy,Pfz) coordinates
        % (x;y;z)_BM = invPf * (x;y;z)
        invPf = arrayfun(@(ii)(inv(Pf(:,:,ii))), 1:size(Pf,3),'UniformOutput',false);
        invPf = cat(3,invPf{:});
        % invPf: 3(Pfx,Pfy,Pfz) x 3(x,y,z) x 3(Pix,Piy,Piz) x nMode x nHkl
        
        % Convert Sab from (x,y,z) coordinates to (xBM,yBM,zBM) coordinates
        % invPf * Sab * invPf'
        invPfM = repmat(permute(invPf,[1 2 4 3]),[1 1 nMode 1]);
        Mab    = mmat(mmat(invPfM,Sab),permute(invPfM,[2 1 3 4]));

        invPf = repmat(permute(invPf,[1 2 4 5 3]),[1 1 3 nMode 1]);
        
        % Pab : 3(Pfx,Pfy,Pfz) x 3(Pix,Piy,Piz) x nMode x nHkl
        Pab = permute(sum(invPf.*repmat(permute(Pab,[5 1 2 3 4]),[3 1 1 1 1]),2),[1 3 4 5 2]);
        
        
        % Convert Sab from (x,y,z) coordinates to (xBM,yBM,zBM) coordinates
        % convert (x,y,z) --> (Pfx,Pfy,Pfz)
        %Mab = permute(sum(invPf.*repmat(permute(Sab,[5 1 2 3 4]),[3 1 1 1 1]),2),[1 3 4 5 2]);
        % convert (x,y,z) --> (Pix,Piy,Piz)
        %Mab = permute(sum(invPf.*repmat(permute(Mab,[5 2 1 3 4]),[3 1 1 1 1]),2),[3 1 4 5 2]);
        
       
        % save all neutron cross sections in spectra
        spectra.intP{ii}  = reshape(intP,[3,nMode,sHkl(2:end)]);
        spectra.Pab{ii}   = reshape(Pab,[3,3,nMode,sHkl(2:end)]);
        spectra.Mab{ii}   = reshape(Mab,[3,3,nMode,sHkl(2:end)]);
    end
    % save unpolarised neutron cross section
    spectra.Sperp{ii} = reshape(Sperp,[nMode,sHkl(2:end)]);
end

% remove cell if there are no twins
if nTwin == 1
    spectra.intP  = spectra.intP{1};
    spectra.Pab   = spectra.Pab{1};
    spectra.Mab   = spectra.Mab{1};
    spectra.Sab   = spectra.Sab{1};
    spectra.Sperp = spectra.Sperp{1};
end

spectra.param.n   = n;
spectra.param.uv  = param.uv;
spectra.param.pol = param.pol;

end
