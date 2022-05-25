function addaniso(obj, matrixIdx, varargin)
% assigns anisotropy to magnetic sites
% 
% ### Syntax
% 
% `addaniso(obj, matrixIdx, {atomTypeIdx}, {atomIdx})`
% 
% ### Description
% 
% `addaniso(obj, matrixIdx, {atomTypeIdx}, {atomIdx})` assigns the
% $[3\times 3]$ matrix selected by `matrixIdx` (using either the matrix
% label or matrix index) to the magnetic sites selected by `atomTypeIdx`
% that can contain a name of an atom or its atom index (see [spinw.atom]).
% If `atomTypeIdx` is not defined, anisotropy will be assigned to all
% magnetic atoms.
% 
% ### Examples
% 
% To show the effect of a fourfold axis on anisotropy, we add $A_1$
% easy-axis anisotropy to atoms at $(1/4,1/4,1/2)$ and plot the result. The
% 3D plot shows anistropy using ellipsoid around the magnetic atoms.
%
% ```
% >>cryst = spinw
% >>cryst.genlattice('lat_const',[4 4 3],'spgr','P 4')
% >>cryst.addatom('r',[1/4 1/4 1/2],'S',1)
% >>cryst.addmatrix('label','A1','value',diag([-0.1 0 0]))
% >>cryst.gencoupling
% >>cryst.addaniso('A1')
% >>plot(cryst)
% >>snapnow
% ```
% 
% ### Input arguments
%
% `obj`
% : [spinw] object.
%
% ### Name-Value Pair Arguments
% 
% `matrixIdx`
% : Either an integer, that selects the matrix according to
%   `obj.matrix.mat(:,:,matrixIdx)`, or a string identical to one
%   of the previously defined matrix labels, stored in
%   `obj.matrix.label`.
% 
% `atomTypeIdx`
% : String or cell of strings that select magnetic atoms by
%   their label. Also can be a vector that contains integers, the index of
%   the magnetic atoms in `obj.unit_cell`, with all symmetry equivalent
%   atoms. Maximum value is $n_{atom}$, if undefined anisotropy is assigned to
%   all magnetic atoms. Optional.
%
% `atomIdx`
% : A vector that contains indices selecting some of the
%   symmetry equivalent atoms. Maximum value is the number of symmetry
%   equivalent atoms generated corresponding to `atomTypeIdx` site. If
%   crystal symmetry is not 0, `atomIdx` is not allowed, since the
%   anisotropy matrix for equivalent atoms will be calculated using the
%   symmetry operators of the space group. Optional.
% 
% ### Output Arguments
% 
% The function adds extra entries in the `obj.single_ion.aniso` field of the
% obj [spinw] object.
% 
% ### See Also
% 
% [spinw], [spinw.single_ion], [spinw.addcoupling], [spinw.addg] and [spinw.addmatrix]
%

mAtom    = obj.matom;
nMagAtom = size(mAtom.idx,2);

if nMagAtom==0
    error('spinw:addaniso:NoMagAtom','There is no magnetic atom in the unit cell with S>0!');
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
    error('spinw:addaniso:WrongCouplingTypeIdx','Input matrix does not exists!');
end

if nargin > 3
    atomIdx = varargin{2};
    if size(obj.lattice.sym, 3) > 1
        error('spinw:addaniso:SymmetryProblem','atomIdx is not allowed when crystal symmetry is higher than P1!');
    end

end

addField = obj.single_ion;

if nargin > 2
    atomTypeIdx = varargin{1};
    
    if length(obj.single_ion.aniso) ~= nMagAtom
        addField.aniso = zeros(1,nMagAtom);
    end
    
    % select atoms by label
    if ischar(atomTypeIdx)
        atomTypeIdx = {atomTypeIdx};
    end
    
    if iscell(atomTypeIdx)
        % loop over all atom labels
        isSelectedAtom = zeros(1,obj.natom);
        for ii = 1:numel(atomTypeIdx)
            % atom index
            %sAIdx = find(strcmp(obj.unit_cell.label,atomTypeIdx{ii}));
            findRes = strfind(obj.unit_cell.label,atomTypeIdx{ii});
            sAIdx = find(~cellfun(@(C)isempty(C),findRes));
            
            if numel(sAIdx) == 0
                error('spinw:addaniso:WrongString','Atom label doesn''t exist!')
            %elseif numel(sAIdx) > 1
            %    error('spinw:addaniso:WrongString','Multiple atoms have the same label!')
            end
            isSelectedAtom(sAIdx) = true;
        end
        atomTypeIdx = find(isSelectedAtom);
    end
    
    for ii = 1:length(atomTypeIdx)
        aTemp = addField.aniso(mAtom.idx == atomTypeIdx(ii));
        if nargin > 3
            aTemp(atomIdx) = matrixIdx;
        else
            aTemp = matrixIdx;
        end
        addField.aniso(mAtom.idx == atomTypeIdx(ii)) = aTemp;
    end
else
    addField.aniso = zeros(1,nMagAtom)+matrixIdx;
end

addField.aniso = int32(addField.aniso);
obj.single_ion = addField;

end