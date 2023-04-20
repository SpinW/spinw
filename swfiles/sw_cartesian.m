function [vyOut, vzOut, vxOut] = sw_cartesian(n)
% creates a right handed Cartesian coordinate system
%
% ### Syntax
%
% `[vy, vz, vx] = sw_cartesian(n)`
%
% `V = sw_cartesian(n)`
%
% ### Description
%
% `[vy, vz, vx] = sw_cartesian(n)` creates an $(x,y,z)$ right handed
% Cartesian coordinate system with $v_x$, $v_y$ and $v_z$ defining the
% basis vectors.
%
% `V = sw_cartesian(n)` the generated basis vectors are stored in the `V`
% matrix: `V = [vx vy vz]` as column vectors.
%
% ### Input Arguments
%
% `n`
% : Either a 3 element row/column vector or a $[3\times 3]$ matrix with
%   columns defining 3 vectors.
%
% ### Output Arguments
%
% `vy,vz,vx`
% : Vectors defining the right handed coordinate system. They are
%           either column of row vectors depending on the shape of the
%           input `n`.
%

if nargin == 0
    swhelp sw_cartesian
    return
end

% Shape of original vector.
if numel(n) == 3
    nShape = size(n);
    n = n(:);
    
    z = [0; 0;-1];
    y = [0;-1; 0];
    
    c = cross(n,z);
    if isa(n,'sym')
        opt = any(sw_always(c));
    else
        opt = any(c);
    end
    
    if opt
        vy = c;
    else
        vy = cross(n,y);
    end
    vz = cross(n,vy);
    
elseif all(size(n) == [3 3])
    if det(n) == 0
        error('sw_cartesian:WrongInput','The input vectors are not linearly independent!')
    end
    
    nShape = [3 1];
    vz = cross(n(:,1),n(:,2));
    vy = cross(vz,n(:,1));
    n = n(:,1);
    
else
    error('sw_cartesian:WrongInput','Wrong size of n!')
end


if nargout == 1
    % return a matrix
    vyOut = [n/norm(n) vy/norm(vy) vz/norm(vz)];
else
    % conserve the shape of the input vector.
    vyOut = reshape(vy/norm(vy),nShape);
    vzOut = reshape(vz/norm(vz),nShape);
    vxOut = reshape(n/norm(n),nShape);
end

end