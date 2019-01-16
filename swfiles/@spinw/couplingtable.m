function bonds = couplingtable(obj, varargin)
% creates tabulated list of all bonds as stored
%
% bonds = COUPLINGTABLE(obj,{bondIdx})
%
% Input:
%
% obj       spinw class object. 
% bondIdx   List of bond indices, by default all bonds will be output.
%           Optional. If a bond index is mutiplied by -1, the table output
%           is a matlab built in table type, works only for Matlab R2013b
%           or later versions.
%
% Output:
%
% bonds is a struct type data that contains the following fields:
%   table   Matrix, where every column defines a bond. The rows are the
%           following: (dl_x, dl_y, dl_z, atom1, atom2, idx, mat_idx1,
%           mat_idx2, mat_idx3). Where (dl_x, dl_y, dl_z) defines the
%           translation vector between the origin of the unit cells of the
%           two interacting atom (if they are in the same unit cell, all
%           three components are zero) from atom1 to atom2. atom1 and atom2
%           are the indices of the atoms in the obj.matom list. idx is the
%           index of the bond, where equivalent bonds have identical
%           indices, typically index is increasing with bond length. The
%           last 3 rows (mat_idx) contains pointers to matrices if they
%           are defined, otherwise zeros.
%   bondv   Additional information for every bond defined in the .table
%           field. The first three rows define the vector pointing from
%           atom1 to atom2 in lattice units. The last row define the bond
%           length in Angstrom.
%   matrix  Contains the coupling matrix for every bond, dimensions are
%           [3 3 nCoupling].
%
% Example:
%
% ...
% crystal.gencoupling
% bonds = crystal.couplingtable(-[1 2 3]).table
%
% This will list only the 1st, 2nd and 3rd neighbour bonds in a formatted
% table.
%
% See also SPINW.MATOM, SPINW.INTMATRIX, SPINW.ADDCOUPLING, SPINW.GENCOUPLING.
%

% no bonds are defined
if isempty(obj.coupling.idx)
    bonds.table = zeros(9,0);
    bonds.bondv = zeros(4,0);
    bonds.matrix = zeros(3,3,0);
    return
end

if nargin < 2
    bondIdx = 1:obj.coupling.idx(end);
else
    bondIdx = varargin{1};
    bondIdx = bondIdx(:)';
end

% create the table
table0 = [obj.coupling.dl;obj.coupling.atom1;obj.coupling.atom2;obj.coupling.idx;obj.coupling.mat_idx];
table0 = table0(:,ismember(obj.coupling.idx,int32(bondIdx)));

% create the bond vectors
mAtom = obj.matom.r;
r1 = mAtom(:,table0(4,:));
r2 = mAtom(:,table0(5,:));
bonds.bondv = r2+double(table0(1:3,:))-r1;
% bond distance in Angstrom
bonds.bondv(4,:) = sqrt(sum((obj.basisvector*bonds.bondv).^2,1));

% create the matrices
Jall = obj.intmatrix('extend',false).all;

nCoupling = size(table0,2);

% assign bond index to the coupling
bonds.matrix = zeros(3,3,nCoupling);
if obj.symbolic
    bonds.matrix = sym(bonds.matrix);
end

bondT = table0(1:5,:);
bondT = bondT(:)';

for ii = 1:size(Jall,2)
    idx = strfind(bondT,Jall(1:5,ii)');
    for jj = 1:numel(idx)
        [idxI, idxJ] = ind2sub([5 nCoupling],idx(jj));
        if idxI == 1
            bonds.matrix(:,:,idxJ) = bonds.matrix(:,:,idxJ) + reshape(Jall(6:14,ii),[3 3]);
        end
    end
end

bonds.label = {'dl_x' 'dl_y' 'dl_z' 'atom1' 'atom2' 'idx' 'mat_idx1' 'mat_idx2' 'mat_idx3'};
bonds.table = table0;

end