function [ Shape ] = AllSegmentsD2( Shape )
%AllSegmentsD2 Computes the D2 descriptor of each segment in a segmented
% shape structure.

% % Draw a random seed for later use:
% seed = randi(100000, 1);

nseg = length(Shape.Segments);

for i=1:nseg
    
    V = Shape.Segments{i}.Vertices;
    F = Shape.Segments{i}.Faces;
    
%     V = V * Shape.Segments{i}.BB.rot;
%     n = size(V, 1);
%     s = 1 ./ Shape.Segments{i}.BB.scale;
%     V = V .* repmat(s', n, 1);
    
    mesh.faces = F;
    mesh.vertices = V;

    % Compute both D1 and D2 histograms:
    Shape.Segments{i}.D1 = calc_D1(mesh);
    Shape.Segments{i}.D2 = calc_D2(mesh);
    
end;

end

