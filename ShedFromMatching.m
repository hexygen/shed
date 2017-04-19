function [ shed, costs ] = ShedFromMatching( Shape1, Shape2, Matching, W )
% SHEDFROMMATCHING Computes SHED (shape edit distance) for the given shapes
% according to the given matching and weights.
% 
% 
% Shape1, Shape2 = Input shapes to match. Shapes should be initialized with
%                  a bounding box and descriptors (obtained using
%                  AllSegmentsBoundingBoxes and AllSegmentsDescriptors
%                  functions)
%
% Matching = Binary matrix which provides the matching between shape parts.
%            Computed by MatchShapes.
%
% W = Weights structure that can be learned from an example set and 
%     contains the following parameters. If no weights are given the
%     default weights are used.
%       wGeometry = cost of changing the geometry of a part
%       wScale = cost of changing the scale of a part
%       wPosition = cost of changing the position or connectivity of a part
%       wDuplicate = cost of duplicating a part
%
% Outputs:
% shed = a number that reflects the shape edit distance between the shapes.
%
% costs = a structure that contains the following costs. This can be used
%         to quickly change the weights and recompute SHED without recomputing the
%         costs.
%   Geometry cost: computed for each pair in the matching based on D1/D2 descriptors.
%                  A complete change in shape geometry should cost "1 unit"
%                  Weighted by the average volume of the part in the two shapes.
%
%   Scale cost: computed for each pair of pairs based on their scale / volume.
%               Doubling (or halving) the size of a single part should cost "1 unit"
%               Weighted by the average volume of the four parts in the two pairs.
%       
%   Position cost: computed for each pair of pairs based on their adjacency.
%                  Switching the positions of two parts should cost "1 unit"
%                  Weighted by the average volume of the four parts in the two pairs.
%
%   Duplicate cost: computed for each duplicated part by counting how many 
%                   matches each part have.
%                   Duplicating a single part should cost "1 unit" for each
%                   additional copy.
%                   Weighted by the total volume of the duplicated instances.
%
%%% If you use this code, please cite the following paper:
%  
%  SHED: Shape Edit Distance for Fine-grained Shape Similarity 
%  Yanir Kleiman, Oliver van Kaick, Olga Sorkine-Hornung, Daniel Cohen-Or 
%  SIGGRAPH ASIA 2015
%
%%% Copyright (c) 2015 Yanir Kleiman <yanirk@gmail.com>


if (nargin < 4)
    % Default weights. These weights gave good results for the sets of
    % shapes that appear in the SHED paper.
    W.wGeometry = 0.4795;
    W.wScale = 0.1258;
    W.wPosition = 0.0034;
    W.wDuplicate = 0.3914;
end;

tic;

alpha = 0.5;    % Controls the ratio between D1 and D2 costs
n = size(Matching, 1);
m = size(Matching, 2);

% Find the indices of matches:
[row, col] = find(Matching);
N = length(row);

%% Step 1: Computing Geometry costs - shape similarity between each two parts:
GeometryCosts = zeros(N, 1);
MatchVolume = zeros(N, 1);

for i=1:N
    a = row(i);
    b = col(i);
    Sa = Shape1.Segments{a};
    Sb = Shape2.Segments{b};

    costD2 = chi2_dist(Sa.D2, Sb.D2);
    costD1 = chi2_dist(Sa.D1, Sb.D1);
    
    GeometryCosts(i) = ((1 - alpha) .* costD2 + alpha .* costD1);
    
    MatchVolume(i) = (Sa.BB.volume + Sb.BB.volume) / 2;
end;

totalMatchVolume = sum(MatchVolume);
% Normalize volume based so the sum of all match volumes is 1:
MatchVolume = MatchVolume / totalMatchVolume;

%% Step 2: Computing Scale costs, by measuring the difference in scales in pairs of pairs.

ScaleCosts = zeros(N);
PositionCosts = zeros(N);

for i=1:N-1
    for j=i+1:N
        a = row(i);
        b = col(i);
        c = row(j);
        d = col(j);
        
        % Getting the volume of each part from the structure:
        Sa = Shape1.Segments{a}.BB.volume;
        Sb = Shape2.Segments{b}.BB.volume;
        Sc = Shape1.Segments{c}.BB.volume;
        Sd = Shape2.Segments{d}.BB.volume;
        
        % Those ratios are 1 if the parts are similar in volume:
        if (Sa < Sb)
            % Pairs scale change from a to b (Sab is always less than 1):
            Sab = Sa / Sb;
            Scd = Sc / Sd;
        else
            % Pairs scale change from b to a (Sab is always less than 1):
            Sab = Sb / Sa;
            Scd = Sd / Sc;
        end;
        
        % This value is 1 if the two pairs have similar CHANGE of volume.
        % It is high when the two pairs have very different change of
        % volume. And it is 2 if the change in one pair is exactly twice
        % the change in the other pair.
        Smax = max(Sab/Scd, Scd/Sab);
        
        % Make sure the cost is 1 when one pair have exactly twice the change 
        % in scale the other pair have. The range of values is now 0 to Inf.
        costScale = Smax - 1;
        
        ScaleCosts(i, j) = costScale;
        ScaleCosts(j, i) = ScaleCosts(i, j);
    end;
end;

PairsVolume = MatchVolume * MatchVolume';

%% Step 3: Compute Position costs by comparing the adjacency of each pair of matches.
% If a part has duplicated matches, check the difference in adjacency 
% between the closest pair of duplicated parts.

% First compute the position cost based on adjacency for all pairs:
for i=1:N-1
    for j=i+1:N
        a = row(i);
        b = col(i);
        c = row(j);
        d = col(j);

        if (a ~= c && b ~= d)
            adj_ac = Shape1.Adj(a, c);
            adj_bd = Shape2.Adj(b, d);
            % Cost is zero for pairs with the same adjacency, higher for
            % pairs with different adjacency, and highest for pairs where
            % parts are adjacent in one shape and far away in the other:
            adj_val = (max(adj_ac, adj_bd)) / (min(adj_ac, adj_bd)) - 1;

            PositionCosts(i, j) = adj_val;
            PositionCosts(j, i) = adj_val;
        end;

    end;
end;


% Then fix pairs which are duplicates such that the error is based on the
% closest instance:
% For each row with duplicated parts in shape 2:
for i=1:n
    if (sum(Matching(i, :)) > 1)

        % Find indices that belong to that row:
        ind = find(row == i);
        minCost = min(PositionCosts(ind, :));
        
        PositionCosts(ind, :) = repmat(minCost, length(ind), 1);
        % The cost is no longer symmetric: each duplicated part is compared with
        % its real neighbors and they on the other hand can be compared to
        % the closest instance of the duplicated part.
    end;
end;
    
% For each column with duplicated parts in shape 1:
for j=1:m
    if (sum(Matching(:, j)) > 1)
        
        % Find indices that belong to that column:
        ind = find(col == j);
        minCost = min(PositionCosts(ind, :));
        
        PositionCosts(ind, :) = repmat(minCost, length(ind), 1);
    end;
end;


%% Step 4: Computing Duplication costs, by aggregating the volume of all duplicated parts:
% When there are no duplicated parts the total match volume is the exactly
% the average of the total volumes of the two shapes.

totDuplicate = (totalMatchVolume - Shape1.TotalVolume / 2 - Shape2.TotalVolume / 2) / totalMatchVolume;
if (totDuplicate < 0)
    totDuplicate = 0;
end;

%% Step 5: Aggregate all the costs into one cost function using the given weight parameters:

% Geometry costs are normalized by the relative volume of a match:
totGeometry = sum(GeometryCosts .* MatchVolume);
% Scale and Position costs are normalized by the combined volume of the 
% pair of matches:
totScale = sum(sum(ScaleCosts .* PairsVolume));
totPosition = sum(sum(PositionCosts .* PairsVolume));


costs.GeometryCosts = GeometryCosts;
costs.ScaleCosts = ScaleCosts;
costs.PositionCosts = PositionCosts;
costs.AddedCosts = zeros(N, 1);
% Keeping the MatchVolume and PairsVolume in the structure, mostly for
% debugging:
costs.MatchVolume = MatchVolume;
costs.PairsVolume = PairsVolume;

costs.totGeometry = totGeometry;
costs.totScale = totScale;
costs.totPosition = totPosition;
costs.totDuplicate = totDuplicate;

% Computing the final SHED value:
shed = W.wGeometry * totGeometry + W.wScale * totScale + W.wPosition * totPosition + W.wDuplicate * totDuplicate;




end

