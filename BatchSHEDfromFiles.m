function [ shed, Matchings, Shapes ] = BatchSHEDfromFiles( shapes_id )
%BATCHPMD Calculates pmd for several shapes in a batch, and outputs all of
% the distances, matching, etc.

path_start = 'C:\Yanir\Documents\Projects\SHED\data\Lamps\BigSet\conv\L';
% path_start = 'C:\Yanir\Documents\Projects\SHED\data\candle\All\conv\';
% path_start = 'c:\Yanir\Documents\Projects\SHED\data\Vases\Set1\';
% path_start = 'c:\Yanir\Documents\Projects\SHED\data\Articulated\conv\';
% path_start = 'c:\Yanir\Documents\Projects\SHED\data\Planes\GoodSet2\P';
% off_end = '.off';
% seg_end = '.seg';
off_end = '-graphcut.off';
seg_end = '-graphcut.seg';

n = length(shapes_id);

shed = zeros(n);
Matchings = cell(n);
Shapes = cell(n, 1);

% Read shapes:
tic;

for i=1:n
    off = [path_start num2str(shapes_id(i)) off_end];
    seg = [path_start num2str(shapes_id(i)) seg_end];

    % Parse segmented shape files:
    Shape = ReadSegmentedShape(off, seg);
    % Compute the bounding box of each shape:
    Shape = AllSegmentsBoundingBoxes(Shape);
    Shape = AllSegmentsD2(Shape);
    
    Shapes{i} = Shape;

    display(['Read shape ' num2str(shapes_id(i)) '.']);
end;

t = toc;
display(['Read files and analyze parts in ' num2str(t) ' seconds.']);

% ticID = tic;
% % Compute PMD:
% for i=1:n
%     for j=i+1:n
%         [curr_shed, curr_matching] = ShapesSHED(Shapes{i}, Shapes{j});
% 
%         shed(i, j) = curr_shed;
%         shed(j, i) = curr_shed;
%         Matchings{i, j} = curr_matching;
%         Matchings{j, i} = curr_matching;
%         
%         display(['Computed SHED(' num2str(i) ', ' num2str(j) ').']);
%         
%     end;
% end;

% t = toc(ticID);
% display(['Computed all PMDs in ' num2str(t) ' seconds.']);

end

