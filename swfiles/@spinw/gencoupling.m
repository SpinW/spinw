function gencoupling(obj, varargin)
% generates bond list
%
% ### Syntax
%
% `gencoupling(obj,Name,Value)`
%
% ### Description
%
% `gencoupling(obj,Name,Value)` generates all bonds up to a certain length
% between magnetic atoms. It also groups bonds based either on crystal
% symmetry (is space group is not $P0$) or bond length (with `tolDist`
% tolerance) is space group is not defined. Sorting bonds based on length
% can be forced by setting the `forceNoSym` parameter to true. To check
% whether a space group is defined call the [spinw.symmetry] function.
%
% {{warning This function has to be used after the crystal structure is defined.
%   The [spinw.addcoupling] function will only work afterwards. }}
%
% ### Examples
%
% A triangular lattice is generated using `spinw.gencoupling` and
% the [spinw.table] function lists the 1st, 2nd and 3rd neighbor bonds:
%
% ```
% >>cryst = spinw
% >>cryst.genlattice('lat_const',[3 3 5],'angled',[90 90 120])
% >>cryst.addatom('r',[0 0 0],'S',1)
% >>cryst.gencoupling
% >>cryst.table('bond',1:3)>>
% ```
%
% ### Input Arguments
%
% `obj`
% : [spinw] object.
%
% ### Name-Value Pair Arguments
%
% `'forceNoSym'`
% : If `true`, equivalent bonds are always generated based on
%   bond length with `tolDist` length tolerance and effectively reducing
%   the bond symmetry to `P0`. If `false` symmetry operators will be used
%   if they are given ([spinw.symmetry] returns `true`).
%
% `'maxDistance'`
% : Maximum bond length that will be stored in the
%   [spinw.coupling] property in units of \\ang. Default value is 8.
%
% `'maxSym'`
% : Maximum bond length until the symmetry equivalent bonds are
%   generated. It is usefull if long bonds have to be generated for the
%   dipolar interaction, but the symmetry analysis of them is not
%   necessary. Default value is equal to `maxDistance`.
%
% `'tolDist'`
% : Tolerance of distance, within two bonds are considered
%   equivalent, default value is $10^{-3}$\\ang. Only used, when no
%   space group is defined.
%
% `'tolMaxDist'`
% : Tolerance added to maxDistance to ensure bonds between same atom in
%   neighbouring unit cells are included when maxDistance is equal to a
%   lattice parameter.
%
% `'dMin'`
% : Minimum bond length, below which an error is triggered.
%   Default value is 0.5 \\ang.
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
% The [spinw.coupling] field of `obj` will store the new bond list, while
% overwriting previous bond list. This will also remove any previous
% assignment of exchange matrices to bonds.
%
% ### See Also
%
% [spinw] \| [spinw.symmetry] \| [spinw.nosym]
%

% is there any symmetry operator?
isSym = size(obj.lattice.sym,3) > 0;

inpForm.fname  = {'forceNoSym' 'maxDistance' 'tolMaxDist' 'tolDist' 'dMin' 'maxSym' 'fid'};
inpForm.defval = {false         8             1e-5  1e-3      0.5   []       -1   };
inpForm.size   = {[1 1]        [1 1]         [1 1] [1 1]     [1 1]  [1 1]   [1 1] };
inpForm.soft   = {false        false          false false    false  true    false };


param = sw_readparam(inpForm, varargin{:});
pref = swpref;

if any([param.maxDistance, param.tolMaxDist, param.tolDist, param.dMin, param.maxSym] < 0)
    error('spinw:gencoupling:NegativeDistances','All distances should be positive.');
end

tolMaxDist   = param.tolMaxDist;
tolD  = param.tolDist;

if tolD > param.maxDistance
    error('spinw:gencoupling:TolDist', ...
          ['Tolerance on symmetrically equivalent bond lengths ' ...
           'larger than max distance.']);
end

% to avoid some problem with symmetry if maxDistance equal to a lattice
% constant
if tolMaxDist > param.maxDistance
    error('spinw:gencoupling:TolMaxDist', ...
        'Tolerance on max distance larger than max distance.');
end
param.maxDistance = param.maxDistance + tolMaxDist;

if param.maxDistance < param.dMin
    error('spinw:gencoupling:MaxDLessThanMinD', ...
        'maxDistance is smaller then dMin parameter');
end


if isempty(param.maxSym)
    param.maxSym = param.maxDistance;
end

if param.fid == -1
    fid = pref.fid;
else
    fid = param.fid;
end

% calculate the height of the paralellepiped of a unit cell
hMax1 = 1./sqrt(sum(inv(obj.basisvector').^2,2));

% calculate the distance of the [1 1 1] point from the origin
hMax2 = abs(sum(obj.basisvector,2));

% calculate the closest point to the origin
hMax = min([hMax1 hMax2],[],2);

% gives the number of extra unit cells along all 3 axes
% that are necessary to cover the minimum bond distance
nC = ceil(param.maxDistance./hMax);

% force no symmetry operator mode
if param.forceNoSym
    isSym = false;
end

if fid ~= 0
    fprintf0(fid,['Creating the bond list (maxDistance = %g ' obj.unit.label{1}...
        ', nCell = %dx%dx%d)...\n'],param.maxDistance-tolMaxDist,nC);
end

% save the sym/nosym method into obj
obj.sym = isSym;

% Generate all atomic positions in the unit cell.
mAtom = obj.matom;

% Number of magnetic atoms in the unit cell.
nMagAtom = size(mAtom.r,2);

if nMagAtom == 0
    error('spinw:gencoupling:NoMagAtom','There is no magnetic atom (S>0) in the unit cell!');
end

% Use half 'cube' around the center unit cell and remove the identical
% bonds using unique() to simplify the calculation

% number of unit cells in the half 'cube'
cDim = [nC(1)+1 2*nC(2)+1 2*nC(3)+1];
nHalfCube = prod(cDim);

% generate all cell translations
cTr{1} = 0:nC(1);
cTr{2} = -nC(2):nC(2);
cTr{3} = -nC(3):nC(3);

% generate cell origin positions in l.u.
[cTr{1}, cTr{2}, cTr{3}] = ndgrid(cTr{:});
% cell origin translations: Na x Nb x Nc x 1 x 1 x3
cTr = cat(6,cTr{:});
% remove unit cells that would produce duplicate bonds
% mark them with NaN (enough to do along a-axis values)
cTr(1,:,(1:end)<=nC(3),1,1,1)       = nan;
cTr(1,(1:end)<=nC(2),nC(3)+1,1,1,1) = nan;
% positions of mgnetic atoms in the half 'cube' in l.u.
% Na x Nb x Nc x 1 x nMagAtom x 3
% atom2
aPos = bsxfun(@plus,permute(mAtom.r,[3:6 2 1]),cTr);
% generate all distances from the atoms in the (0,0,0) cell in l.u.
% r(atom2) - r(atom1)
% Na x Nb x Nc x nMagAtom x nMagAtom x 3
dR  = bsxfun(@minus,aPos,permute(mAtom.r,[3:5 2 6 1]));
% mark duplicate bonds within the (0,0,0) cell with nan
R0 = permute(dR(1,nC(2)+1,nC(3)+1,:,:,:),[4:6 1:3]);
dR(1,nC(2)+1,nC(3)+1,:,:,:) = bsxfun(@times,tril(nan(nMagAtom))+1,R0);
% calculate the absolute value of the distances in Angstrom
dRA = sqrt(sum(mmat(permute(obj.basisvector,[3:7 1 2]),dR,[6 7]).^2,6));
% reshape the numbers into a column list of bonds
% 3 x Na x Nb x Nc x nMagAtom x nMagAtom
dl    = permute(repmat(cTr,[1 1 1 nMagAtom nMagAtom 1]),[6 1:5]);
atom1 = repmat(permute(1:nMagAtom,[1 3:5 2 6]),[1 cDim 1 nMagAtom]);
atom2 = repmat(permute(1:nMagAtom,[1 3:5 6 2]),[1 cDim nMagAtom 1]);
dRA   = permute(dRA,[6 1:5]);
% store everything in a single matrix
cMat  = [dl(:,:);atom1(1,:);atom2(1,:);dRA(1,:)];
% remove nan bonds
cMat = cMat(:,~isnan(cMat(1,:)));
% sort according to increasing distance
% cMat = sortrows(cMat',6)'; too slow
[~,cIdx] = sort(cMat(6,:));
cMat     = cMat(:,cIdx);
% check cutoff for bond lengths
if cMat(6,1) < param.dMin
    error('spinw:gencoupling:AtomPos',['Some atoms are too close (Dmin=' ...
        num2str(cMat(6,1)) '<' num2str(param.dMin) '), check your crystal structure!']);
elseif cMat(6,1) > param.maxDistance
    error('spinw:gencoupling:maxDistance', ...
        'There are no bonds with distance < maxDistance');
end

cMat = cMat(:,cMat(6,:)<=param.maxDistance);
% index the bonds
cMat(7,:) = cumsum([1 diff(cMat(6,:))>tolD]);

dRA  = cMat(6,:);
% keep the bond length
cMat = cMat([1:5 7],:);
% cMat rows: [la, lb, lc, atom1, atom2, idx]

% symmetry equivalent couplings
nSym = int32(0);
if isSym
    idxSym = dRA<=param.maxSym;
    if ~any(idxSym)
        warning('spinw:gencoupling:maxSym',['No bonds with distances < '...
            'maxSym = ', num2str(param.maxSym), ' - all bonds will be ' ...
            'reduced to P0 symmetry.']);
    else
        % store the final sorted couplings in nMat
        nMat = zeros(6,0);
        ii  = 1;
        idx = 1;
        % maximum bond index for symmetry operations
        maxidxSym = max(cMat(6,idxSym));

        while ii <= maxidxSym
            % select columns from sorM with a certain idx value
            sortMs = cMat(:,cMat(6,:) == ii);
            while size(sortMs,2)>0
                [genC, unC] = swsym.bond(mAtom.r, obj.basisvector, sortMs(:,1), obj.lattice.sym, tolD);
                genCAll = [genC [-genC(1:3,:);genC([5 4],:)]];
                % remove from sortMs the identical couplings
                iNew = isnew(genCAll(1:5,:),sortMs(1:5,:), tolD);
                sortMs(:,~iNew) = [];
                % remove identical couplings from the symmetry generated
                % list
                genC(:,~unC) = [];
                if sum(~iNew) ~= sum(unC)
                    error('spinw:gencoupling:SymProblem', ...
                        ['Symmetry error! ii=%d, idx=%d. Try to '...
                        'change ''tol'' parameter.'], ii, idx)
                end
                % move the non-unique (not new) couplings (symmetry equivalent ones)
                nMat = [nMat [genC;ones(1,size(genC,2))*idx]]; %#ok<AGROW>
                idx  = idx + 1;
            end
            ii = ii + 1;
        end

        % include the increase of bond index in case bonds are splitted due to
        % symmetry inequivalency
        cMat = cMat(:,cMat(6,:)>maxidxSym);
        cMat(6,:) = cMat(6,:) + nMat(6,end)-maxidxSym;
        cMat = [nMat cMat];

        % save the value of maximum bond index that is generated by symmetry
        nSym = int32(nMat(6,end));
    end
end

% default anisotropy and g-tensor values
aniso = int32(zeros(1,nMagAtom));
g     = int32(zeros(1,nMagAtom));

if fid ~= 0
    fprintf0(fid,'...%d bonds are retained out of %d generated!\n',size(cMat,2),nHalfCube*nMagAtom^2);
end

% save output structure
coupling.dl      = int32(cMat(1:3,:));
coupling.atom1   = int32(cMat(4,:));
coupling.atom2   = int32(cMat(5,:));
coupling.mat_idx = int32(zeros(3,size(cMat,2)));
coupling.idx     = int32(cMat(6,:));
coupling.type    = int32(zeros(3,size(coupling.idx,2)));
coupling.sym     = int32(zeros(3,size(coupling.idx,2)));
coupling.rdip    = obj.coupling.rdip;
coupling.nsym    = nSym;

obj.coupling         = coupling;
obj.single_ion.aniso = aniso;
obj.single_ion.g     = g;

spinw.validate(obj);

end

function [isnew, symIdx] = isnew(A,B, tol)
% selects the new vectors from B. Dimensions of A and B have to be [3 nA]
% and [3 nB] respectively.
%
% A vector in B is considered new, if d(mod(vA-vB,1))<tol.
%

nA = size(A,2);
nB = size(B,2);

notequal = sum(abs(repmat(permute(A,[2 3 1]),[1 nB 1]) - repmat(permute(B,[3 2 1]),[nA 1 1])).^2,3) > tol;

isnew = all(notequal,1);

idx = 1:nB;

symIdx = arrayfun(@(idx)find(~notequal(:,idx),1,'first'),idx(~isnew));

end