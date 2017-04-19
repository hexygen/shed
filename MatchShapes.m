function [ matching, Aff_out, Shape1, Shape2 ] = MatchShapes( Shape1, Shape2 )
%%% MatchShapes Computes a matching between the given shapes parts.
% Each shape construct is given as following:
%   - Shape.Segments contains a cell array of segments, where each segement
%     contains Nx3 matrix of vertices in that segment.
%   - Shape.Adj contains an MxM extended adjacency matrix of the segments.
%     This matrix indicates for each pair of segments (i, j) how far they 
%     are from each other.
%   - Shape.BB contains the bounding box of the entire shape for volume
%     calculations.
%   
% The Shape objects should inlcude the bounding box and D2 computations of
% each part (i.e. this should be called after AllSegmentsDescriptors).
%
%
%%% If you use this code, please cite the following paper:
%  
%  SHED: Shape Edit Distance for Fine-grained Shape Similarity 
%  Yanir Kleiman, Oliver van Kaick, Olga Sorkine-Hornung, Daniel Cohen-Or 
%  SIGGRAPH ASIA 2015
%
%%% Copyright (c) 2015 Yanir Kleiman <yanirk@gmail.com>


tic;

% alpha controls the relative weight of D1 and D2 in the geometry cost.
% Value in SHED paper: alpha = 0.5 (equal weights)
alpha = 0.5;
% beta controls the weight of second-order terms compared to first-order
% terms.
% Value in SHED paper: beta = 0.3
beta = 0.3;
% gamma controls the relative weight of adjacency and scale in the
% second-order term.
% Value in SHED paper: gamma = 0.5 (equal weights)
gamma = 0.5;
% sigma is used to transform costs into affinities using an inverse
% exponent.
% Value in SHED paper: sigma = 0.5
sigma = 0.5;

n = length(Shape1.Segments);
m = length(Shape2.Segments);


%% Step 1: Computing shape similarity between each two parts:
PairsCostD2 = zeros(n, m);
PairsCostD1 = zeros(n, m);

for i=1:n
    for j=1:m
        Si = Shape1.Segments{i};
        Sj = Shape2.Segments{j};

        % Alternative ways to compare the histograms are commented out:
%         costD2 = earthMoverDist(Si.D2, Sj.D2);
%         costD2 = norm(Si.D2 - Sj.D2).^2;
        costD2 = chi2_dist(Si.D2, Sj.D2);
        PairsCostD2(i, j) = costD2;

%         costD1 = earthMoverDist(Si.D1, Sj.D1);
%         costD1 = norm(Si.D1 - Sj.D1).^2;
        costD1 = chi2_dist(Si.D1, Sj.D1);
        PairsCostD1(i, j) = costD1;
        
    end;
end;

PairsCost = ((1 - alpha) .* PairsCostD2 + alpha .* PairsCostD1);

%% Step 2: Computing first and second order affinities between matches
% First order affinities are the geometric similarity of two parts in a
% match.
% Second order affinities are computed according to the adjacency of parts
% in both matches and their scaling similarity.

nm = n*m;
Costs = zeros(nm, nm);  % Costs matrix
Aff = zeros(nm, nm);    % Affinity matrix
matsize = [n, m];


for ind1 = 1:nm
    % Getting parts index:
    [a, b] = ind2sub(matsize, ind1);

    % First order costs = geometric similarity between parts:
    cost = PairsCost(a, b);
    Costs(ind1, ind1) = cost;

    % Affinity is derived directly from cost using a gaussian kernel:
    Aff(ind1, ind1) = exp(-cost/sigma);

    % Second order costs are computed between two matches, ind1 and ind2
    % or (a, b) and (c, d):
    for ind2 = ind1+1:m*n
        [c, d] = ind2sub(matsize, ind2);


        adj_ac = Shape1.Adj(a, c);
        adj_bd = Shape2.Adj(b, d);

        % The pairs are considered compatible if they have a similar
        % adjacency and similar change of scale between them:
        sa = Shape1.Segments{a}.BB.volume;
        sb = Shape2.Segments{b}.BB.volume;
        sc = Shape1.Segments{c}.BB.volume;
        sd = Shape2.Segments{d}.BB.volume;

        sac = sa / sc;
        sbd = sb / sd;

        if (adj_ac == 0 || adj_bd == 0)
            aff_val = 0;
        else
            adj_val = (max(adj_ac, adj_bd)) / (min(adj_ac, adj_bd)) - 1;
            % adj_val = 0 when the adjacencies match, and high when the
            % adjacencies do not match. Also higher when the parts are
            % close: for example 5/4 < 3/2.
            
            scale_val = (max(sac, sbd)) / (min(sac, sbd)) - 1;

            % Weighted average of adjacency and scaling factor:
            val = gamma * adj_val + (1-gamma) * scale_val;
            exp_val = exp(-val);

            % Relative weight of second-order terms compared to first-order
            % terms:
            aff_val = beta * exp_val;
        end;
        
        % Assign second-order terms to affinity matrix:
        Aff(ind1, ind2) = aff_val;
        Aff(ind2, ind1) = aff_val;
    end;
end;

% Output original affinity before the matching algorithm which updates the
% matrix:
Aff_out = Aff;

%% Step 3: Iterative, adaptive spectral matching algorithm
% In each step, the eigenvector is computed, the highest value is picked
% and incosistencies are added to the objective matrix for the next
% eigenvector computation.
matching = zeros(matsize);
match_ind = zeros(nm, 1);

done = 0;
ind_max = 0;
exclude = [];

% Iterative, adaptive spectral matching algorithm:
while (~done)
    % Get principal eigenvector of the affinity matrix:
    [x, ~] = eigs(Aff, 1);
    % Make sure the vector is positive and avoid rounding errors:
    x = abs(x);
    
    % If a match was selected it shouldn't be selected again:
    x(match_ind == 1) = 0;
    x(exclude) = 0;
    
    % Find maximum value of eigenvector - most likely pair to match:
    [x_max, ind_max] = max(x);

    % Find sub-indices of maximum value:
    [a, b] = ind2sub(matsize, ind_max);

    % summing the rows to get an n*1 vector:
    sn = sum(matching, 2);
    % summing the columns to get an m*1 vector:
    sm = sum(matching, 1);
    
    % If every part is matched with at least one other part, we are done:
    if ((min(sn) > 0 && min(sm) > 0) || x_max == 0)
        done = 1;
    else
    
        % This match is only possible if one of the parts is not yet
        % matched:
        if ((sn(a) == 0) || (sm(b) == 0))

            matching(a, b) = 1;
            match_ind(ind_max) = 1;

            % Update affinity matrix:
            ind1 = sub2ind(matsize, a * ones(1, m), 1:m);
            ind2 = sub2ind(matsize, 1:n, b * ones(1, n));

            % If one pair contains part a and the other contains part b they
            % are not compatible:
            Aff(ind1, ind2) = 0;
            Aff(ind2, ind1) = 0;

            % The affinity of the selected match becomes 1 to strengthen
            % any compatible matches:
            Aff(ind_max, ind_max) = 1;

            % The affinity between all previously selected pairs is 1 since
            % they are already selected:
            Aff(logical(match_ind), logical(match_ind)) = 1;
        else
            % The match is not valid given previous matches, add it to the 
            % excluded matches list:
            exclude = [exclude ind_max];
        end;
        
    end;
end;

%% Remove incorrect duplications such as many-to-many relations in the matching graph:
% This should not happen in practice, but just to be safe...

% summing the rows to get an n*1 vector:
sn = sum(matching, 2);
% summing the columns to get an m*1 vector:
sm = sum(matching, 1);

for a=1:n
    for b=1:m
        
        if (matching(a, b) && sn(a) > 1 && sm(b) > 1)
            % Many-to-many relation found!
            % Arbitrarily remove the first one we encounter:
            matching(a, b) = 0;
            
            
            % summing the rows to get an n*1 vector:
            sn = sum(matching, 2);
            % summing the columns to get an m*1 vector:
            sm = sum(matching, 1);            
        end;
    end;
end;



t = toc;

display(['Matching computed in ' num2str(t) ' seconds.']);

end

