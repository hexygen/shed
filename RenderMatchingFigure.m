function [ ] = RenderMatchingFigure( Shapes, Matchings, filename, shapes_list, output_dir, rot_step )
% RenderMatchingFigures renders a nice figure of the matching between 
% shapes.
% Each segment in the source shape is assigned a unique color, which is 
% used to render both the original segment and any matching segment in the 
% target shape. If several segments in the source shape match the same 
% segment in the target shape, the segment will be rednered several times
% with random perturbations so that it will appear in several different 
% colors.
%
% Input:
% Shapes = the collection of shapes (can include more shapes than the ones
%          that are compared)
% Matchings = a cell array where cell {i, j} is a matching matrix between
%             shape i and shape j
% filename = Prefix of the file name to save for each figure. 
%            The figure of the source shape will be saved as
%            "filename_Segs".
%            The figures of the target shapes will be saved as 
%            "filename_Matching_id.png" where id is the index of the shape 
%            in the collection.
% shapes_list = a list of indices of shapes. The first index is the source
%               shape and the rest are the target shapes to compare with.
% output_dir = a directory name in which to save the figures.
%
% rot = rotation of each figure around the Y axis, for better viewing angle.
%
%
% Output:
% n files will be saved into the designated folder.
%
%%% If you use this code, please cite the following paper:
%  
%  SHED: Shape Edit Distance for Fine-grained Shape Similarity 
%  Yanir Kleiman, Oliver van Kaick, Olga Sorkine-Hornung, Daniel Cohen-Or 
%  SIGGRAPH ASIA 2015
%
%%% Copyright (c) 2015 Yanir Kleiman <yanirk@gmail.com>

src_id = shapes_list(1);
target_ids = shapes_list(2:end);
k = length(target_ids); % number of target shapes

% This is the gap between segments, relative to the distance from the shape
% center. Change to 0 to have a render without any gaps.
gap = 0.05;

if (nargin < 6)
    % Show only one rotation:
    rot_step = 0;
    rot_num = 1;
else
    % Show many rotations:
    rot_num = floor(2 * pi / rot_step);
end;

% Add a fixed rotation to avoid straight-on viewing angle:
fix_rot = pi / 20;

% For vases, candles, lamps:
elevation = -pi/10;

% For planes:
% elevation = -pi/5;

% This color map is used in order to color segments.
% Colors: [Yellow, Blue, Red, Turquoise, Green, Orange, Purple, Dark Blue, 
%          Dark Green, Dark Red, Light Gray, Dark Gray, Pink]
Colors = [0.9 0.9 0.5;
          0.15 0.6 0.8;
          1 0.4 0.4;
          0.4 0.8 0.7;
          0.5 0.8 0.45;
          1 0.5 0.33;
          0.55 0.33 0.7;
          0.35 0.45 1;
          0.2 0.8 0.5;
          0.8 0.4 0.2;
          0.75 0.75 0.75;
          0.5 0.5 0.5;
          0.95 0.65 0.8];


%% Render source shape:
seg = Shapes{src_id}.Segments;
s = length(seg);    % number of segments

shape_center = Shapes{src_id}.BB.center;

% Produce K rotations of each image:
for r=1:rot_num

    figure('Name', 'Source Shape');
    axis off;


    % For each segment:
    for i=1:s
        mesh = struct();
        mesh.vertices = seg{i}.Vertices;
        mesh.faces = fliplr(seg{i}.Faces);

        % Rotate mesh for better looking figures:
        mesh = mesh_rotate(mesh, [0 1 0], rot_step * r + fix_rot); %0.678 for lamps and vases, 0.12 for candles
        mesh = mesh_rotate(mesh, [1 0 0], elevation); 

        % Add relative distance from shape center:
        seg_center = mean(mesh.vertices);
        delta = (seg_center - shape_center) * gap;
        delta = repmat(delta, length(mesh.vertices), 1);
        mesh.vertices = mesh.vertices + delta;


        % Plot segment:
        p = patch(mesh);

        % Color segment:
        set(p, 'FaceColor', Colors(i, :), 'EdgeColor', 'none', 'BackFaceLighting', 'lit');
        set(p, 'AmbientStrength', 0.8);


    end;

    % After all patches were added to the figure, apply render settings to 
    % make the figure nice:
    MakeFigureNice();

    % Save current figure:
    saveas(gcf, [output_dir filesep filename '_Segs_' num2str(r) '.png'], 'png');

    % Close current figure:
    close(gcf);
end;

%% Render target shapes:
for j=1:k
    
    
    seg = Shapes{target_ids(j)}.Segments;
    shape_center = Shapes{target_ids(j)}.BB.center;
    st = length(seg);
    
    matching = Matchings{src_id, target_ids(j)};

    % Produce K rotations of each image:
    for r=1:rot_num

        figure('Name', ['Target Shape #' num2str(j) ', r=' num2str(r)]);
        axis off;
        
        % Go over the original segments and find their matches in the current shape:
        for i=1:s
            mult = 1;
            % For each matching segment:
            for ii=1:st
                if (matching(i, ii))

                    % Take patch data from target shape (ii):
                    mesh.vertices = seg{ii}.Vertices;
                    mesh.faces = fliplr(seg{ii}.Faces);

                    % Add relative distance from shape center:
                    seg_center = mean(mesh.vertices);
                    delta = (seg_center - shape_center) * gap;
                    delta = repmat(delta, length(mesh.vertices), 1);
                    mesh.vertices = mesh.vertices + delta;
                    
                    % Add random perturbation to the segments' vertices:
                    mesh.vertices = mesh.vertices + (rand(size(mesh.vertices)) - 0.5) / 1000;

                    % Rotate mesh for better looking figures:
                    mesh = mesh_rotate(mesh, [0 1 0], rot_step*r + fix_rot); %0.678 for lamps and vases, 0.12 for candles
                    mesh = mesh_rotate(mesh, [1 0 0], elevation); 

                    % Plot segment:
                    p = patch(mesh);

                    % Set color according to source shape (i):
                    set(p, 'FaceColor', Colors(i, :) * mult, 'EdgeColor', 'none', 'BackFaceLighting', 'lit');
                    set(p, 'AmbientStrength', 0.8);
                    
                    mult = mult * 0.94;
                end;
            end;
        end;

        % After all patches were added to the figure, apply render settings to 
        % make the figure nice:
        MakeFigureNice();
        
        
        % Save current figure:
        saveas(gcf, [output_dir filesep filename '_Matching_'  num2str(target_ids(j)) '_' num2str(r) '.png'], 'png');

        % Close current figure:
        close(gcf);
        
    end;
    
end;



% This function makes a figure look nice in terms of lighting etc.
function [ ] = MakeFigureNice( )
    
% Set aspect ratio
daspect([1 1 1]);

material dull;

% Set lightning
lighting gouraud;

%%% Three point lighting (sort of)!
light('Position',[0.6 0.6 1],'Style','infinite', 'Color', [0.7 0.7 0.7]);
light('Position',[-1 0.3 -1],'Style','infinite', 'Color', [0.25 0.25 0.25]);
light('Position',[0 -1 1],'Style','infinite', 'Color', [0.3 0.3 0.3]);
light('Position',[-1 1 0.2],'Style','infinite', 'Color', [0.45 0.45 0.45]);
light('Position',[-0.5 0.3 1],'Style','infinite', 'Color', [0.7 0.7 0.7]);


% BLUISH Ambient:
set(gca, 'AmbientLightColor', [0.45 0.57 0.65]);


% CLAY Ambient:
% set(gca, 'AmbientLightColor', [0.5 0.36 0.3]);
% set(p, 'AmbientStrength', 1);



% Set view
view(0, 90);
    
end


end

