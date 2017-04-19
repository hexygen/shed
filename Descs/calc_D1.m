function d1 = calc_D1(mesh)
%%% This function calculates D1 descriptor of a mesh.
% See: Shape Distributions, Osada et al.
%
% The shape is sampled according to the respective area of each face.
%
% Written by Noa Fish, 2013

n = 1024;
b = 64;

faces = mesh.faces;
vertices = mesh.vertices;

% Normalize vertices:
vn = size(vertices, 1);
center = mean(vertices);
vertices = vertices - repmat(center, vn, 1);

% Compute areas:
areas = computeArea(mesh);

cum_areas = cumsum(areas);
cum_areas = cum_areas - cum_areas(1);
cum_areas = cum_areas / cum_areas(end);

dists = zeros(n, 1);
p = zeros(n, 3);

r = rand(n, 3);

index = arrayfun(@(z) sum(z >= cum_areas), r(:, 1));
vind = faces(index, :);
a1 = vertices(vind(:, 1), :);
b1 = vertices(vind(:, 2), :);
c1 = vertices(vind(:, 3), :);
r1 = r(:, 2);
r2 = r(:, 3);
o1 = [1 1 1];
p = (1-sqrt(r1))*o1.*a1 + (sqrt(r1).*(1-r2))*o1.*b1 + (sqrt(r1).*r2)*o1.*c1;

dists = sqrt(sum(p.^2, 2));


d1 = hist(dists, b) / n;


end