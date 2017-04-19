function [ shed, costs ] = BatchShedFromMatching( Shapes, Matchings, Weights )
% BATCHSHEDFROMMATCHING Given a collection of shapes and the matchings
% between them, computes pairwise SHED for the entire collection.
%
% Shapes = the collection of shapes (cell array)
%
% Matchings = a cell array where each cell {i, j} contains the matching
%             matrix of shape i to shape j
%
% Weights = Weights structure that can be learned from an example set and 
%           contains the following parameters. If no weights are given the
%           default weights are used.
%       wGeometry = cost of changing the geometry of a part
%       wScale = cost of changing the scale of a part
%       wPosition = cost of changing the position or connectivity of a part
%       wDuplicate = cost of duplicating a part
%
%
% Output:
% shed = a matrix of the SHED value of each pair of shapes.
%
% costs = a collection of costs, where each cost is a matrix of the final
%         cost for each pair of shapes. This is useful to know which costs
%         contributed most to the final SHED value for each pair of shapes.
%         Contains costGeometry, costScale, costPosition and costDuplicate.
%
%%% If you use this code, please cite the following paper:
%  
%  SHED: Shape Edit Distance for Fine-grained Shape Similarity 
%  Yanir Kleiman, Oliver van Kaick, Olga Sorkine-Hornung, Daniel Cohen-Or 
%  SIGGRAPH ASIA 2015
%
%%% Copyright (c) 2015 Yanir Kleiman <yanirk@gmail.com>


n = length(Shapes);

shed = zeros(n);
% Most code here takes care of outputing the costs...
costs.costGeometry = zeros(n);
costs.costScale = zeros(n);
costs.costPosition = zeros(n);
costs.costDuplicate = zeros(n);

for i=1:n-1;
    for j=i+1:n;
        
        if (nargin > 2)
            % Use given weights:
            [curr_shed, cost] = ShedFromMatching(Shapes{i}, Shapes{j}, Matchings{i, j}, Weights);
        else
            % Use default weights:
            [curr_shed, cost] = ShedFromMatching(Shapes{i}, Shapes{j}, Matchings{i, j});
        end;
        
        shed(i, j) = curr_shed;
        
        costs.costGeometry(i, j) = cost.totGeometry;
        costs.costScale(i, j) = cost.totScale;
        costs.costPosition(i, j) = cost.totPosition;
        costs.costDuplicate(i, j) = cost.totDuplicate;
    end;
end;

% Make matrices symmetric (not sure this part is necessary):
shed = shed + shed';

costs.costGeometry = costs.costGeometry + costs.costGeometry';
costs.costScale = costs.costScale + costs.costScale';
costs.costPosition = costs.costPosition + costs.costPosition';
costs.costDuplicate = costs.costDuplicate + costs.costDuplicate';

end

