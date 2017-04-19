%
% This function rotates a mesh by a certain angle on a certain axis
%
% function N = mesh_rotate(M, axis, angle)
%
% Input -
%   - M: triangle mesh: M.vertices(i, :) represents the 3D coordinates
%   of vertex 'i', while M.faces(i, :) contains the indices of the three
%   vertices that compose face 'i'
%    - axis: rotation axis: a 3D vector
%    - angle: rotation angle in radians: a scalar
%
% Output -
%   - N: new rotated triangle mesh, with the same structure as M
%
function N = mesh_rotate(M, axis, angle)
%
% Copyright (c) 2008 Oliver van Kaick <ovankaic@cs.sfu.ca>
%

% Create output mesh
N = struct();
N.vertices = M.vertices;
N.faces = M.faces;

% Rotate vertex positions
R = build_rotation_matrix(axis, angle);
for i = 1:size(M.vertices, 1)
    N.vertices(i,:) = M.vertices(i,:)*R;
end

% Copy colors, if defined
if (isfield(M, 'FaceVertexCData'))
    N.FaceVertexCData = M.FaceVertexCData;
end


%
% This function builds a generalized rotation matrix
%
% Input -
%   - axis: rotation axis: a 3D vector
%   - angle: rotation angle: a scalar
%
% Output -
%   - R: rotation matrix: a 3x3 matrix
%
function R = build_rotation_matrix(axis, angle)

% Init matrix
R = zeros(3,3);

% Normalize axis vector
axis = axis / norm(axis, 2);
x = axis(1);
y = axis(2);
z = axis(3);

% Compute sin and cos
cosa = cos(angle);
sina = sin(angle);

% Create matrix
% http://en.wikipedia.org/wiki/Rotation_matrix
R(1,1) = cosa + (1 - cosa)*x^2;
R(1,2) = (1 - cosa)*x*y - sina*z;
R(1,3) = (1 - cosa)*x*z + sina*y;

R(2,1) = (1 - cosa)*y*x + sina*z;
R(2,2) = cosa + (1 - cosa)*y^2;
R(2,3) = (1 - cosa)*y*z - sina*x;

R(3,1) = (1 - cosa)*z*x - sina*y;
R(3,2) = (1 - cosa)*z*y + sina*x;
R(3,3) = cosa + (1 - cosa)*z^2;
