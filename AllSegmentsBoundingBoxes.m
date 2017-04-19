function [ Shape ] = AllSegmentsBoundingBoxes( Shape )
% AllSegmentsBoundingBoxes Finds the bounding box of each segment in a 
% segmented shape structure.
%
%%% If you use this code, please cite the following paper:
%  
%  SHED: Shape Edit Distance for Fine-grained Shape Similarity 
%  Yanir Kleiman, Oliver van Kaick, Olga Sorkine-Hornung, Daniel Cohen-Or 
%  SIGGRAPH ASIA 2015
%
%%% Copyright (c) 2015 Yanir Kleiman <yanirk@gmail.com>

nseg = length(Shape.Segments);

Shape.TotalVolume = 0;
Shape.TotalRootVolume = 0;

for i=1:nseg
    % Calculates the bounding box of the current segment:
    BB = CalcBoundingBox(Shape.Segments{i}.Vertices);
    Shape.Segments{i}.BB = BB;
    % Total volume of the shape is the sum of the volume of all segments:
    Shape.TotalVolume = Shape.TotalVolume + BB.volume;
end;

Shape.AverageVolume = Shape.TotalVolume / nseg;

end

