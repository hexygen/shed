% ------------------------------------------------------------------------
%  SHED: Shape Edit Distance for Fine-grained Shape Similarity 
%  Yanir Kleiman, Oliver van Kaick, Olga Sorkine-Hornung, Daniel Cohen-Or 
%  SIGGRAPH ASIA 2015
% ------------------------------------------------------------------------
%
% This code computes Shape Edit Distance (SHED) 
% for a collection of segmented shapes.
% It is well documented, so feel free to take a look inside.
%
% First, add the folder 'Descs' and its subfolder into the path.
% The starting point is the function ShedFromList.
% I attached 7 vases as an example to get you going.
% To test SHED on those vases, run the following command:

[shed, Matchings, Shapes, shape_files] = ShedFromList('list.txt');
 
% You can use the output values to further experiment with the system,
% for example by changing the weights, providing a manual matching, etc.
% 
% The most important part of the code is:
% MatchShapes.m: The adaptive spectral matching algorithm described in the paper.
% ShedFromMatching.m: Computes SHED from loaded segmented shapes with a given matching.
%
% To render the results, use the following code:

RenderMatchingFigure(Shapes, Matchings, 'test', 1:7, 'testdata');

% It renders and saves a figure for each shape. The segments are color
% coded to match the first shape in the list (in this case, shape #1).
% The figures are saved in the folder 'data' (the last parameter).
% Look inside the function to see more options.
%
% To segment shapes using approximately convex segmentation, see here:
% http://www.cs.tau.ac.il/~noafish/wcseg/
%
%%% Copyright (c) 2015 Yanir Kleiman <yanirk@gmail.com>


