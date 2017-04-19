function [  ] = ShowMatch( Si, Sj, name )
% ShowMatch Shows a match between segment i and segment j in one figure.
%
% Input:
% Si, Sj = Segments from the Shape construct, for example:
%          ShowMatch(Shape.Segments{i}, Shape2.Segments{j}, figure_name);
% name = Name of the output figure.
%
%%% If you use this code, please cite the following paper:
%  
%  SHED: Shape Edit Distance for Fine-grained Shape Similarity 
%  Yanir Kleiman, Oliver van Kaick, Olga Sorkine-Hornung, Daniel Cohen-Or 
%  SIGGRAPH ASIA 2015
%
%%% Copyright (c) 2015 Yanir Kleiman <yanirk@gmail.com>

figure('Name', name);

if (~isempty(Si))

    v = Si.Vertices;
    subplot(1, 2, 1);
    scatter3(v(:, 1), v(:, 2), v(:, 3), 30, 'filled'); 
    axis equal;
end;
    
if (~isempty(Sj))

    v = Sj.Vertices;
    subplot(1, 2, 2);
    scatter3(v(:, 1), v(:, 2), v(:, 3), 30, 'filled'); 
    axis equal;
end;

end

