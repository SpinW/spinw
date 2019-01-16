Aniso = diag([-1 0 0]);

%% you define 2 vectors (the easy axis of Aniso (100) above and the easy axis direction you want in your system (111)):

v1 = [1 0 0];
v2 = [1 1 1]/sqrt(3);

%Find the axis around which v1 can be rotated to v2:

ax = cross(v1,v2);

%Find the rotation angle:

phi = atan2(norm(cross(v1,v2)),dot(v1,v2));

%Create a rotation matrix:

R = sw_rotmat(ax,phi);

% Rotate the Aniso matrix:

A = R*Aniso*R;

%This A matrix will define an easy axis anisotropy along the (111) direction with the size of 1 meV. To double check that it is right, we calculate the eigenvalues:

[V, Aniso2] = eig(A);
