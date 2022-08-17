function atomList = atom(obj)
% generates symmetry equivalent atomic positions
% 
% ### Syntax
% 
% `atomList = atom(obj)`
% 
% ### Description
% 
% `atomList = atom(obj)` generates all atomic positions using the symmetry
% operators stored in `obj.lattice.sym` and the symmetry inequivalent
% atomic positions in `obj.unit_cell.r`. If no symmetry is defined (denoted
% $P0$ symmetry) or the symmetry is $P1$ the function returns simply the
% positions stored in `obj.unit_cell.r`.
% 
% ### Examples
% 
% Here we create a new space group, that contains all the translations of
% the FCC lattice. Then create a crystal with an atom at `[0 0 0]` position.
% The `cryst.atom` lists all 4 symmetry equivalent positions generated using
% the FCC symmetry operators:
%
% ```
% >>cryst = spinw
% >>opStr = 'x+1/2,y+1/2,z;x+1/2,y,z+1/2;x,y+1/2,z+1/2';
% >>cryst.genlattice('lat_const',[8 8 8],'sym',opStr,'label','FCC')
% >>cryst.addatom('r',[0 0 0],'label','Atom1')
% >>atomList = cryst.atom
% >>atomList.r>>
% ```
%
% ### Input Arguments
% 
% `obj`
% : [spinw] object.
% 
% ### Output Arguments
% 
% `atomList` is a structure with the following fields:
% * `r`     Positions of the atoms in lattice units stored in matrix with
%           dimensions of $[3\times n_{atom}]$. 
% * `idx`   Indices of the atoms in the [spinw.unit_cell] field stored in a
%           matrix with dimensions of $[1\times n_{atom}]$.
% * `mag`   Vector of logical variables, `true` if the spin of the atom is
%           non-zero, dimensions are $[1\times n_{atom}]$.
% 
% ### See Also
% 
% [spinw] \| [spinw.matom] \| [swsym.add] \| [spinw.genlattice] \| [spinw.addatom]
%

% Generate all symmetry equivalent atoms
[atomList.r, atomList.idx] = swsym.position(obj.lattice.sym,obj.unit_cell.r);

if obj.symbolic
    atomList.mag = ~sw_always(obj.unit_cell.S(atomList.idx)==0);
else
    atomList.mag = obj.unit_cell.S(atomList.idx)>0;
end

end