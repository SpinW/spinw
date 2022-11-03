function optm = optmagsteep(obj, varargin)
% quench optimization of magnetic structure
% 
% ### Syntax
% 
% `optm = optmagsteep(obj,Name,Value)`
% 
% ### Description
% 
% `optm = optmagsteep(obj,Name,Value)` determines the lowest energy
% magnetic configuration within a given magnetic supercell and previously
% fixed propagation (and normal) vector (see [spinw.optmagk]). It
% iteratively rotates each spin towards the local magnetic field thus
% achieving local energy minimum. Albeit not guaranteed this method often
% finds the global energy minimum. The methods works best for small
% magnetic cells and non-frustrated structures. Its execution is roughly
% equivalent to a thermal quenching from the paramagnetic state.
% 
% ### Input Arguments
% 
% `obj`
% : [spinw] object.
% 
% ### Name-Value Pair Arguments
% 
% `'nRun'`
% : Number of iterations, default value is 100 (it is usually enough). Each
%   spin will be quenched `nRun` times or until convergence is reached.
% 
% `'boundary'`
% : Boundary conditions of the magnetic cell, string with allowed values:
%   * `'free'`  Free, interactions between extedned unit cells are
%               omitted.
%   * `'per'`   Periodic, interactions between extended unit cells
%               are retained.
%
%   Default value is `{'per' 'per' 'per'}`.
% 
% `'nExt'`
% : The size of the magnetic cell in number of crystal unit cells.
%   Default value is taken from `obj.mag_str.nExt`.
% 
% `'fSub'`
% : Function that defines non-interacting sublattices for parallelization.
%   It has the following header:
%       `cGraph = fSub(conn,nExt)`, where `cGraph` is a row vector with
%       $n_{magExt}$ number of elements,
%   `conn` is a matrix with dimensions of $[2\times n_{conn}]$ size matrix and $n_{ext}$ is equal to
%   the `nExt` parameter. Default value is `@sw_fsub`.
% 
% `'subLat'`
% : Vector that assigns all magnetic moments into non-interacting
%   sublattices, contains a single index $(1,2,3...)$ for every magnetic
%   moment in a row vector with $n_{magExt}$ number of elements. If
%   undefined, the function defined in `fSub` will be used to partition the
%   lattice.
% 
% `'random'`
% : If `true` random initial spin orientations will be used (paramagnet),
%   if initial spin configuration is undefined (`obj.mag_str.F` is empty)
%   the initial configuration will be always random. Default value is
%   `false`.
% 
% `'TolX'`
% : Minimum change of the magnetic moment necessary to reach convergence.
% 
% `'saveAll'`
% : Save moment directions for every loop, default value is `false`.
% 
% `'Hmin'`
% : Minimum field value on the spin that moves the spin. If the
%   molecular field absolute value is below this, the spin won't be
%   turned. Default is 0.
% 
% `'plot'`
% : If true, the magnetic structure in plotted in real time. Default value
%   is `false`.
% 
% `'pause'`
% : Time in second to pause after every optimization loop to slow down plot
%   movie. Default value is 0.
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
% `optm`
% : Struct type variable with the following fields:
%   * `obj`         spinw object that contains the optimised magnetic structure.
%   * `M`           Magnetic moment directions with dimensions $[3\times n_{magExt}]$, if
%                   `saveAll` parameter is `true`, it contains the magnetic structure
%                   after every loop in a matrix with dimensions $[3\times n{magExt}\times n_{loop}]$.
%   * `dM`          The change of magnetic moment vector averaged over all moments
%                   in the last loop.
%   * `e`           Energy per spin in the optimised structure.
%   * `param`       Input parameters, stored in a struct.
%   * `nRun`        Number of loops executed.
%   * `datestart`   Starting time of the function.
%   * `dateend`     End time of the function.
%   * `title`       Title of the simulation, given in the input.
% 
% ### See Also
% 
% [spinw] \| [spinw.anneal] \| [sw_fsub] \| [sw_fstat]
%

% disable warning in spinw.energy
warnStruct = warning('off','spinw:energy:AnisoFieldIncomm');

% save the time of the beginning of the calculation
if nargout > 0
    optm.datestart  = datestr(now);
end

% get magnetic structure
nExt   = double(obj.mag_str.nExt);
title0 = 'Optimised magnetic structure using the method of steepest descent';

inpForm.fname  = {'nRun' 'epsilon' 'random' 'boundary'          'subLat' 'Hmin' 'pause'};
inpForm.defval = {100     1e-5      false   {'per' 'per' 'per'}  []      0      0      };
inpForm.size   = {[1 1]   [1 1]     [1 1]   [1 3]                [1 -1]  [1 1]  [1 1]  };
inpForm.soft   = {0       0         0       0                    1       false  false  };

inpForm.fname  = [inpForm.fname  {'nExt' 'fSub'   'TolX' 'title' 'saveAll' 'plot' 'fid'}];
inpForm.defval = [inpForm.defval {nExt   @sw_fsub 1e-10  title0  false     false  -1   }];
inpForm.size   = [inpForm.size   {[1 3]  [1 1]    [1 1]  [1 -2]  [1 1]     [1 1]  [1 1]}];
inpForm.soft   = [inpForm.soft   {0      0        0      0        0        false  false}];


param = sw_readparam(inpForm,varargin{:});
pref = swpref;

if prod(param.nExt) == 0
    error('spinw:optmagsteep:WrongInput','''nExt'' has to be larger than 0!');
end

% Text output file
if param.fid == -1
    fid = pref.fid;
else
    fid = param.fid;
end

fprintf0(fid,['Optimising the magnetic structure using local spin '...
    'updates\n(nRun = %d, boundary = (%s,%s,%s))...\n'],param.nRun,param.boundary{:});

% Creates random spin directions if param.random is true.
mag_param = struct;
if param.random || isempty(obj.mag_str.F) || any(param.nExt~=nExt)
    mag_param.mode = 'random';
    mag_param.nExt = param.nExt;
    obj.genmagstr(mag_param);
    % TODO check
    nExt = param.nExt;
end

% get the magnetic structure
magStr = obj.magstr('exact',false);

M = magStr.S;

% Produce the interaction matrices
[SS, SI] = obj.intmatrix;

% add the dipolar interactions to SS.all
SS.all = [SS.all SS.dip];


% express translations in the original unit cell
SS.all(1:3,:) = bsxfun(@times,SS.all(1:3,:),nExt');

% Function options.
nRun    = param.nRun;
nMagExt = size(M,2);

% Spin length for normalization.
S     = sqrt(sum(M.^2,1));

% Modify the interaction matrices according to the boundary conditions.
for ii = 1:3
    if strcmp('free',param.boundary{ii})
        SS.all(:,SS.all(ii,:)~=0) = [];
    end
end

% Spins are not allowed to be coupled to themselves. Remove these couplings
% and give a warning.
idxSelf = SS.all(4,:)==SS.all(5,:);
if any(idxSelf)
    warning('spinw:optmagsteep:SelfCoupling','Some spins are coupled to themselves in the present magnetic cell!');
    SS.all(:,idxSelf) = [];
end

% Calculates the energy of the initial configuration and prepares the
% anisotropy matrix. B is in units of the couplings.
Bloc = permute(mmat(SI.field*obj.unit.muB,SI.g),[2 3 1]);
AA = SI.aniso;

% convert all anisotropy matrix to have a maximum eigenvalue of zero
% anisotropy gan be generated from eigenvalues and eigenvectors: A = V*E*V';
for ii = 1:size(AA,3)
    [AAv,AAe]  = eig(AA(:,:,ii));
    AAe2       = diag(diag(AAe)-max(AAe(:)));
    AA(:,:,ii) = AAv*AAe2*AAv';
end

Ax = squeeze(AA(:,1,:));
Ay = squeeze(AA(:,2,:));
Az = squeeze(AA(:,3,:));

% Checks whether there is any external field
param.isfield = any(Bloc(:));

% Checks whether anisotropy is non-zero.
if any(AA(:))
    param.aniso = true;
else
    param.aniso = false;
end

% Assing moments to sublattices for parallel calculation. There are no
% coupling between moments on the same sublattice, thus Weiss field can be
% calculated parallel. SSc stores the index of the sublattice, size:
% (1,nMagExt)

if isempty(param.subLat)
    % add anisotropies
    %Aidx = (squeeze(sumn(abs(AA),[1 2]))>0)';
    %SSc  = param.fSub([SS.all(4:5,:) [Aidx;Aidx]],param.nExt);
    SSc  = param.fSub(SS.all(4:5,:),param.nExt);
    param.subLat = SSc;
else
    SSc  = param.subLat;
end
nSub = max(SSc);

nNeighG = zeros(nMagExt,1);
for ii = 1:nMagExt
    %nNeighG(ii) = sum((SS.all(4,:) == ii)|(SS.all(5,:) == ii));
    nNeighG(ii) = sum(SS.all(4,:) == ii) + sum(SS.all(5,:) == ii);
end

% Maximum number of neighbours
maxNeighG = max(nNeighG);

% Interaction matrices and neigbor indices
SSiG = zeros(maxNeighG,nMagExt) + (nMagExt+1);
SSJG = zeros(9,maxNeighG,nMagExt);

% indices of the transpose of the J matrix in a column of SS.all
% indices are between 6-14
trIdx = reshape(reshape(1:9,[3 3])',[9 1])+5;

% magnetic ordering wave vector
km = magStr.k;

if any(km) && numel(km)>3
    warning('spinw:optmagsteep:Multik','Multi-k structures cannot be optimized!')
    return
end

% for non-zero km, rotate the exchange matrices that couple spins between
% different unit cell
if any(km)
    % Rotate the coupling matrices that couple spins in different unit cells
    % Si * Jij * Sj' = Si * Jij * R * Sj
    % Sj' = R(km,dl) * Sj
    [~,R] = sw_rot(magStr.n,km*SS.all(1:3,:)*2*pi);
    Jrot  = mmat(reshape(SS.all(6:14,:),3,3,[]),R);
    JJR   = reshape(Jrot,9,[]);
    [~,R] = sw_rot(magStr.n,-km*SS.all(1:3,:)*2*pi);
    Jrot  = mmat(reshape(SS.all(trIdx,:),3,3,[]),R);
    JJTR  = reshape(Jrot,9,[]);
else
    JJR   = SS.all(6:14,:);
    JJTR  = SS.all(trIdx,:);
end

% Indexes for transposing J for exchanged spins in the interaction.
% Default is Si * J * Sj, or Sj * J' * Si has to be used.
%trIdx = reshape(reshape(1:9,[3 3])',[9 1])+5;
for ii = 1:nMagExt
    idx1 = SS.all(4,:) == ii;
    idx2 = (SS.all(5,:) == ii);
    SSiG(1:nNeighG(ii),ii)   = [SS.all(5,idx1)    SS.all(4,idx2)    ]';
    %SSJG(:,1:nNeighG(ii),ii) = [SS.all(6:14,idx1) SS.all(trIdx,idx2)];
    SSJG(:,1:nNeighG(ii),ii) = [JJR(:,idx1) JJTR(:,idx2)];
end

% Store spin indices of each sublattice for speedup.
Sindex = zeros(nSub,nMagExt);
Sindex(nSub*(0:(nMagExt-1))+SSc) = 1;
Sindex      = logical(Sindex);

% Remove uncoupled moments, they should keep their original orientation
fSpin = squeeze(sumn(abs(AA),[1 2]))==0 & nNeighG==0 & sum(abs(Bloc'),2)==0;
Sindex(:,fSpin) = false;

if ~any(Sindex)
    error('spinw:optmagsteep:NoField','There nothing to optimise!');
end

% Speeds up the code by storing every sublattice data in different cells
csSSiG = cell(nSub,1);
csSSJG = cell(nSub,1);
cAx    = cell(nSub,1);
cAy    = cell(nSub,1);
cAz    = cell(nSub,1);
cS     = cell(nSub,1);
cB     = cell(nSub,1);

% store sublattice indices in cell
sSindexF = cell(1,nSub);

for ii = 1:nSub
    sSindex    = Sindex(ii,:);
    cS{ii}    = S(sSindex);
    
    cAx{ii}   = Ax(:,sSindex);
    cAy{ii}   = Ay(:,sSindex);
    cAz{ii}   = Az(:,sSindex);
    
    csSSiG{ii} = reshape(SSiG(:,sSindex),1,[]);
    csSSJG{ii} = reshape(SSJG(:,:,sSindex),3,3,[]);
    
    cB{ii}    = Bloc(:,sSindex);
    sSindexF{ii} = find(Sindex(ii,:));
end

if fid == 1
    sw_timeit(0,1,'Magnetic structure optimization');
end

if nargout == 1
    E   = zeros(1,nRun);
    if param.saveAll
        Msave = zeros(3,nMagExt);
    end
end

% Initial step size is infinite, to make at least 1 cycle.
dM = inf;
% Initial index.
rIdx = 0;

% add extra zero moment as a placeholder
M = [M zeros(3,1)];

% create swplot figure if it doesn't exist
if param.plot
    try
        hFigure = swplot.activefigure;
    catch
        hFigure = obj.plot();
    end
end

while (rIdx < nRun) && (dM>param.TolX)
    Mold = M;
    for jsub = 1:nSub
        % Logical vector, selecting the moments on a given
        % sublattice [1,nMagExt]
        %sSindex = Sindex(jsub,:);
        % F stores the molecular field acting on the moments of
        % the jsub sublattice (exchange+external field).
        % F [3, nElementSub(jsub)]
        sSSJG = csSSJG{jsub};
        F = squeeze(sum(reshape(permute(mmat(sSSJG,permute(M(:,csSSiG{jsub}),[1 3 2])),[1 3 2]),3, maxNeighG,[]),2));
        
        % Adds external magnetic field.
        if param.isfield
            F = F - cB{jsub};
        end
        
        % Adds anisotropy field.
        if param.aniso
            % Select the moment vectors on the sublattice.
            %Ms = M(:,[sSindex false]);
            Ms = M(:,sSindexF{jsub});
            Fa = 2*[sum(Ms.*cAx{jsub},1); sum(Ms.*cAy{jsub},1); sum(Ms.*cAz{jsub},1)];
            F = F + Fa;
        end
        
        % molecular field absolute value
        normF = sqrt(sum(F.^2));
        
        if param.Hmin > 0
            % don't move moments where the |F| is smaller than a minimum
            % value
            nzF     = normF >= param.Hmin;
            SindexS = sSindexF{jsub}(nzF);
            F       = F(:,nzF);
            cSS     = cS{jsub}(nzF);
        else
            SindexS = sSindexF{jsub};
            cSS     = cS{jsub};
        end
        %M(:,[sSindex false]) = -bsxfun(@times,F,cS{jsub}./sqrt(sum(F.^2)));
        if ~isempty(F)
            % turn the moments toward the local field
            M(:,SindexS) = -bsxfun(@times,F,cSS./normF);
        end
        
    end
    
    % Calculates the system energy at the end of the temperature step.
    if nargout > 0 || param.plot
        Mexport = M(:,1:(end-1));
        obj.mag_str.F    = Mexport + 1i*cross(repmat(permute(magStr.n,[2 3 1]),[1 nMagExt 1]),Mexport);
        obj.mag_str.nExt = int32(nExt);
        obj.mag_str.k    = km';
    end
    
    if nargout > 0
        E(rIdx+1) = obj.energy;
        if param.saveAll
            Msave(:,:,rIdx+1) = M(:,1:end-1);
        end
    end
    
    % plot magnetic structure
    if param.plot
        swplot.plotmag('figure',hFigure);
        drawnow;
    end

    % Check stopping condition, give the dM limit.
    dM = sum(sqrt(sum((Mold - M).^2,1)))/nMagExt;
    rIdx = rIdx + 1;

    if param.pause > 0
        % wait a bit
        pause(param.pause);
    end
    
    if fid == 1
        sw_timeit(rIdx/param.nRun*100);
    end

end

if fid == 1
    sw_timeit(100,2);
else
    if fid ~= 0
        fprintf0(fid,'Calculation finished.\n');
    end
end

if rIdx == nRun
    warning('spinw:optmagsteep:NotConverged', ...
            'Convergence was not reached!')
end

% Save optimised magnetic structure into the spinw object.
%obj.genmagstr('mode','helical','S',M(:,1:end-1),'k',km,'n',magStr.n,'nExt',nExt);
M = M(:,1:(end-1));
obj.mag_str.F    = M + 1i*cross(repmat(permute(magStr.n,[2 3 1]),[1 nMagExt 1]),M);
obj.mag_str.k    = km';
obj.mag_str.nExt = int32(nExt);


% Create output structure.
if nargout > 0
    optm.obj      = copy(obj);
    if param.saveAll
        optm.M = Msave;
    else
        optm.M = M;
    end
    optm.dM       = dM;
    optm.e        = E(1:rIdx);
    optm.param    = param;
    optm.nRun     = rIdx;
    optm.dateend  = datestr(now);
    optm.title    = param.title;
end

% restore warnings
warning(warnStruct);

end