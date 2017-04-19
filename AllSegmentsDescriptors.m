function [ Shape ] = AllSegmentsDescriptors( Shape )
% AllSegmentsDescriptors Computes the D1 and D2 descriptors of each segment
% in a segmented shape structure.
%
%%% If you use this code, please cite the following paper:
%  
%  SHED: Shape Edit Distance for Fine-grained Shape Similarity 
%  Yanir Kleiman, Oliver van Kaick, Olga Sorkine-Hornung, Daniel Cohen-Or 
%  SIGGRAPH ASIA 2015
%
%%% Copyright (c) 2015 Yanir Kleiman <yanirk@gmail.com>

nseg = length(Shape.Segments);

for i=1:nseg
    
    V = Shape.Segments{i}.Vertices;
    F = Shape.Segments{i}.Faces;
    
    mesh.faces = F;
    mesh.vertices = V;

    % Compute both D1 and D2 histograms:
    Shape.Segments{i}.D1 = calc_D1(mesh);
    Shape.Segments{i}.D2 = calc_D2(mesh);
    
end;

end

