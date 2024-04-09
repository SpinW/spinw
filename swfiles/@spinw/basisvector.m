function basisVector = basisvector(obj, varargin)
% generates lattice vectors
% 
% ### Syntax
% 
% `basisVec = basisvector(obj, {norm})`
% 
% ### Description
% 
% `basisVec = basisvector(obj, {norm})` returns the lattice vectors of the
% unit cell in a $[3\times 3]$ matrix, with the $a$, $b$ and $c$ vectors
% stored in columns. The vectors are normalized to the lattice parameters
% by default.
% 
% ### Examples
% 
% The `basisVec` matrix can be used to change coordinate system, converting
% between positions expressed in lattice units to positions expressed in
% \\ang, using `r_lu` for lattice unit coordinates and `r_xyz` for
% \\ang units (both stored in a column vector) the conversions are the
% following:
% ```
% r_xyz = basisVec * r_lu
% ```
% or
% ```
% r_lu = inv(basisVec)*r_xyz
% ```
%
% It is also possible to convert between momentum vector in reciprocal
% lattice units (rlu) into \\ang$^{-1}$ units. Assuming that momentum
% vectors are row vectors:
% ```
% Q_xyz =  Q_rlu * 2*pi*inv(basisVec)
% ```
% or
% ```
% Q_rlu = 1/(2*pi)*Q_xyz*basisVect
% ```
% 
% ### Input Arguments
% 
% `obj`
% : [spinw] object.
% 
% `norm`
% : If `true`, the basis vectors are normalized to 1, otherwise the
%   length of the basis vectors are equal to the lattice constants. Default
%   is `false`.
% 
% ### Output Arguments
% 
% `basisVec`
% : Stores the three lattice vectors in columns, dimensions are $[3\times 3]$.
% 
% ### See Also
% 
% [spinw] \| [spinw.abc] \| [spinw.rl]
%

if nargin < 2
    norm = false;
else
    norm = varargin{1};
end

if nargin < 3
    symbolic = false;
else
    symbolic = varargin{2};
end

if symbolic
    ang = obj.lattice.angle;

    % check for special angles (multiples of 30)
    specIdx = abs(round(ang*180/pi/30)-ang*180/pi/30)<1e-8;
    angSym = sym([0 0 0]);
    angSym(specIdx) = sym(round(ang(specIdx)*180/pi))*pi/180;
    
    if any(~specIdx)
        symAng0 = [sym('alpha','positive') sym('beta','positive') sym('gamma','positive')];
        angSym(~specIdx) = symAng0(~specIdx);
    end
    
    alpha = angSym(1);
    beta  = angSym(2);
    gamma = angSym(3);
else
    alpha = obj.lattice.angle(1);
    beta  = obj.lattice.angle(2);
    gamma = obj.lattice.angle(3);
end

v1 = [1          0          0];
v2 = [cos(gamma) sin(gamma) 0];

v3(1) = cos(beta);
v3(2) = sin(beta)*(cos(alpha)-cos(beta)*cos(gamma))/(sin(beta)*sin(gamma));
v3(3) = sqrt(sin(beta)^2-v3(2)^2);

basisVector = [v1' v2' v3'];

% if not normalized, use the lattice constants
if ~norm
    basisVector = basisVector*diag(obj.lattice.lat_const);
end

end
