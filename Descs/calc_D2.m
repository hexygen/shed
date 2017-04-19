function d2 = calc_D2(mesh)
%%% This function calculates D2 descriptor of a mesh.
% See: Shape Distributions, Osada et al.
%
% The shape is sampled according to the respective area of each face.
%
% Written by Noa Fish, 2013

n = 1024;
b = 64;

faces = mesh.faces;
vertices = mesh.vertices;
areas = computeArea(mesh);

cum_areas = cumsum(areas);
cum_areas = cum_areas - cum_areas(1);
cum_areas = cum_areas / cum_areas(end);

dists = zeros(n,1);

for i=1:n
    r = rand(2,1);
    
    index1 = binarysearch(cum_areas,r(1));
    vert1 = faces(index1,:);
    a1 = vertices(vert1(1),:);
    b1 = vertices(vert1(2),:);
    c1 = vertices(vert1(3),:);
    r2 = rand(2,1);
    p1 = (1-sqrt(r2(1)))*a1 + sqrt(r2(1))*(1-r2(2))*b1 + sqrt(r2(1))*r2(2)*c1;
    
    index2 = binarysearch(cum_areas,r(2));
    vert2 = faces(index2,:);
    a2 = vertices(vert2(1),:);
    b2 = vertices(vert2(2),:);
    c2 = vertices(vert2(3),:);
    r2 = rand(2,1);
    p2 = (1-sqrt(r2(1)))*a2 + sqrt(r2(1))*(1-r2(2))*b2 + sqrt(r2(1))*r2(2)*c2;
    
    d = norm(p1 - p2);
    dists(i) = d;
end

d2 = hist(dists, b) / n;


end