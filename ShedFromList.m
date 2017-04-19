function [ shed, Matchings, Shapes, shape_files ] = ShedFromList( list_filename )
% SHEDfromList loads all shapes in the list file, and computes the matching
% and SHED between each pair of shapes.
% Note that this function uses DEFAULT weights. To run with different
% weights use BatchShedFromMatching on the output.
%
% Input:
% list_filename = a list of shape files wihtout extension - extensions for
%                 model is ".off" and extension for segments is ".seg".
%
% Output:
% shed = matrix of SHED values for each pair of shapes in the list.
%
% Matchings = a cell array where cell {i, j} is the matching matrix of
%             shape i to shape j
% Shapes = the collection of loaded segmented shapes with their computed
%          bounding boxes and volumes
% shape_files = the list of file names. Useful to see which shape
%               corresponds to each file without opening the file
%               externally.
%
%%% If you use this code, please cite the following paper:
%  
%  SHED: Shape Edit Distance for Fine-grained Shape Similarity 
%  Yanir Kleiman, Oliver van Kaick, Olga Sorkine-Hornung, Daniel Cohen-Or 
%  SIGGRAPH ASIA 2015
%
%%% Copyright (c) 2015 Yanir Kleiman <yanirk@gmail.com>


off_end = '.off';
seg_end = '.seg';

tic;

% Open list file:
fid = fopen(list_filename);
if (fid == -1)
    disp(['ERROR: could not open list file "' list_filename '"']);
    return;
end;

% Read shape file names:
shape_files = {};
while (~feof(fid))
    filename = fgetl(fid);
    if (~isempty(filename))
        shape_files = [shape_files; filename];
    end;
end;

fclose(fid);

t = toc;
display(['Read file names in ' num2str(t) ' seconds.']);

n = length(shape_files);

shed = zeros(n);
Shapes = cell(n, 1);

% Read shapes:
ticID = tic;

for i=1:n
    tic;
    % Add appropriate extensions:
    off = [shape_files{i} off_end];
    seg = [shape_files{i} seg_end];

    % Parse segmented shape files:
    Shape = ReadSegmentedShape(off, seg);
    % Compute the bounding box for each segment in the shape:
    Shape = AllSegmentsBoundingBoxes(Shape);
    % Compute geometric descriptors for each segment:
    Shape = AllSegmentsDescriptors(Shape);
    
    % Output shape:
    Shapes{i} = Shape;

    t = toc;
    display(['Read shape ' shape_files{i} ' in ' num2str(t) ' seconds.']);
end;

t = toc(ticID);
display(['Read files and analyze parts in ' num2str(t) ' seconds.']);

% Find matching between each pair of shapes:
Matchings = BatchMatch(Shapes);

% Compute SHED for each pair of shapes:
shed = BatchShedFromMatching(Shapes, Matchings);

end

