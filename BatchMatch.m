function [ Matchings ] = BatchMatch( Shapes, ind )
% BatchMatch Gets a list of loaded segmented shapes and calculated the
% matchings between all pairs of shapes using MatchShapes().
%
% Shapes = collection of shapes
% ind = the indices in the collection to compare (optional)
%
%%% If you use this code, please cite the following paper:
%  
%  SHED: Shape Edit Distance for Fine-grained Shape Similarity 
%  Yanir Kleiman, Oliver van Kaick, Olga Sorkine-Hornung, Daniel Cohen-Or 
%  SIGGRAPH ASIA 2015
%
%%% Copyright (c) 2015 Yanir Kleiman <yanirk@gmail.com>

if (nargin < 2)
    % Default index list is whole collection:
    ind = 1:length(Shapes);
end;

n = length(ind);
Matchings = cell(n);

ticID = tic;
for i=1:n
    for j=i:n
        % Find matching between shapes i and j:
        curr_matching = MatchShapes(Shapes{ind(i)}, Shapes{ind(j)});

        Matchings{i, j} = curr_matching;
        Matchings{j, i} = curr_matching';
        
        display(['Computed Matching(' num2str(i) ', ' num2str(j) ')']);
        % Uncomment this to show a figure for each matching:
%         ShowMatchingOneFigure(Shapes{ind(i)}, Shapes{ind(j)}, Matchings{i, j});
    end;
end;

t = toc(ticID);
display(['Computed all matchings in ' num2str(t) ' seconds.']);


end

