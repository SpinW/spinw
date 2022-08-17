function addg(obj, matrixIdx, varargin)
% assigns g-tensor to magnetic atoms
% 
% ### Syntax
% 
% `addg(obj, matrixIdx, {atomTypeIdx}, {atomIdx})`
% 
% ### Description
% 
% `addg(obj, matrixIdx, {atomTypeIdx}, {atomIdx})` assigns the
% $[3\times 3]$ matrix selected by `matrixIdx` (using either the matrix
% label or matrix index) to the magnetic sites selected by `atomTypeIdx`
% that can contain a name of an atom or its atom index (see [spinw.atom]).
% If `atomTypeIdx` is not defined, g-tensor will be assigned to all
% magnetic atoms.
% 
% ### Examples
% 
% The following example will add the $g_1$ diagonal matrix to all magnetic
% atoms as anisotropic g-tensor and show the effect of a fourfold axis:
% 
% ```
% >>cryst = spinw
% >>cryst.genlattice('lat_const',[4 4 3],'sym','P 4')
% >>cryst.addatom('r',[1/4 1/4 1/2],'S',1)
% >>cryst.addmatrix('label','g_1','value',diag([2 1 1]))
% >>cryst.gencoupling
% >>cryst.addg('g_1')
% >>cryst.plot('ionMode','g')
% >>snapnow
% ```
% 
% ### Input Arguments
% 
% `matrixIdx`
% : Either an integer, that selects the matrix
%   `obj.matrix.mat(:,:,matrixIdx)`, or a string identical to one
%   of the previously defined matrix labels, stored in
%   `obj.matrix.label`. Maximum value is $n_{mat}$.
% 
% `atomTypeIdx`
% : String or cell of strings that select magnetic atoms by
%   their label. Also can be a vector that contains integers, the index of
%   the magnetic atoms in [spinw.unit_cell], this will assign the given
%   g-tensor to all symmetry equivalent atoms. Maximum value is $n_{atom}$.
%   If `atomTypeIdx` is not defined, the given g-tensor will be assigned to
%   all magnetic atoms. Optional.
%
% `atomIdx`
% : A vector that contains indices selecting some of the
%   symmetry equivalent atoms. Maximum value is the number of symmetry
%   equivalent atoms corresponding to `atomTypeIdx`. If the crystal
%   symmetry is higher than $P0$, `atomIdx` is not allowed, since the
%   g-tensor for equivalent atoms will be calculated using the symmetry
%   operators of the space group. Optional.
% 
% ### Output Arguments
% 
% The function adds extra entries to the `obj.single_ion.g` matrix.
% 
% ### See Also
% 
% [spinw] \| [spinw.addcoupling] \| [spinw.addaniso] \| [spinw.addmatrix]
%



% assigns g-tensor to magnetic ions
%
% ADDG(obj, matrixIdx, {atomTypeIdx}, {atomIdx})
%
% Input:
%
% matrixIdx     Either an integer, that selects the matrix
%               obj.matrix.mat(:,:,matrixIdx), or a string identical to one
%               of the previously defined matrix labels, stored in
%               obj.matrix.label. Maximum value is nJ.
% atomTypeIdx   String or cell of strings that select magnetic atoms by
%               their label. Also can be a vector that contains integers,
%               the index of the magnetic atoms in obj.unit_cell, with all
%               symmetry equivalent atoms. Maximum value is nAtom, if
%               undefined g-tensor is assigned to all magnetic atoms.
%               Optional.
%  atomIdx      A vector that contains indices selecting some of the
%               symmetry equivalent atoms. Maximum value is the number of
%               symmetry equivalent atoms generated. If crystal symmetry is
%               not 0, atomIdx is not allowed, since the g-tensor for
%               equivalent atoms will be calculated using the symmetry
%               operators of the space group. Optional.
%
% Output:
%
% The function adds extra entries in the 'single_ion.g' field of the obj sw
% object.
%
% Example:
%
% ...
% cryst.addmatrix('label','g1','value',diag([1.8 1.8 2.1]))
% cryst.gencoupling
% cryst.addg('g1')
%
% This will add the 'g1' diagonal matrix to all magnetic atoms as
% anisotropic g-tensor.
%
% See also SPINW, SPINW.ADDCOUPLING, SPINW.ADDANISO, SPINW.ADDMATRIX.
%

mAtom    = obj.matom;
nMagAtom = size(mAtom.idx,2);

if nMagAtom==0
    error('spinw:addg:NoMagAtom','There is no magnetic atom in the unit cell with S>0!');
end

cTemp = -1;
if ischar(matrixIdx)
    for ii = length(obj.matrix.label):-1:1
        if strcmp(obj.matrix.label{ii},matrixIdx)
            cTemp = ii;
        end
    end
    matrixIdx = cTemp;
    
end

if matrixIdx == -1
    error('spinw:addg:WrongCouplingTypeIdx','Input matrix does not exists!');
end

% check matrix has the correct properties for a valid g-tensor
% see Eq. 15.38 in Abragam and Bleaney (1970) EPR of Transition Ions
gtensor = obj.matrix.mat(:, :, matrixIdx);
G = gtensor' * gtensor;
eigvals = eig(G);
if ~issymmetric(G) || ~all(eigvals > 10*eps)
    error('spinw:addg:InvalidgTensor',  ...
        'Check input matrix - g*.g must be symmeteric positive-definite');
end

if nargin > 3
    atomIdx = varargin{2};
    if size(obj.lattice.sym, 3) > 1
        error('spinw:addg:SymmetryProblem','atomIdx is not allowed when crystal symmetry is not P1!');
    end

end

addField = obj.single_ion;

if nargin > 2
    atomTypeIdx = varargin{1};
    
    if length(obj.single_ion.g) ~= nMagAtom
        addField.g = zeros(1, nMagAtom, 'int32');
    end
    
    % select atoms by label
    if ischar(atomTypeIdx)
        atomTypeIdx = {atomTypeIdx};
    end
    
    if iscell(atomTypeIdx)
        % loop over all atom labels
        isSelectedAtom = zeros(1,max(mAtom.idx));
        for ii = 1:numel(atomTypeIdx)
            isSelectedAtom = isSelectedAtom | strcmp(obj.unit_cell.label,atomTypeIdx{ii});
        end
        atomTypeIdx = find(isSelectedAtom);
    end
    
    for ii = 1:length(atomTypeIdx)
        aTemp = addField.g(mAtom.idx == atomTypeIdx(ii));
        if nargin > 3
            aTemp(atomIdx) = matrixIdx;
        else
            aTemp = matrixIdx;
        end
        addField.g(mAtom.idx == atomTypeIdx(ii)) = aTemp;
    end
    
else
    addField.g = zeros(1, nMagAtom, 'int32') + matrixIdx;
end

obj.single_ion = addField;

end