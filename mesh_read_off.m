%
% This function reads an off file and returns a mesh structure
%
% function [M, status] = mesh_read_off(filename)
%
% Input -
%   - filename: name of off file to load
%
% Output -
%   - M: triangle mesh: M.vertices(i, :) represents the 3D coordinates
%   of vertex 'i', while M.faces(i, :) contains the indices of the three
%   vertices that compose face 'i'
%   - status: this variable is 0 if the file was succesfuly opened, or 1
%   otherwise
%
% See also mesh_read
%
function [M, status] = mesh_read_off(filename)
%
% Copyright (c) 2008, 2009 Oliver van Kaick <ovankaic@cs.sfu.ca>
%

    % Open file
    fid = fopen(filename);
    status = 0;
    if fid == -1
        disp(['ERROR: could not open file "' filename '"']);
        M = [];
        status = 1;
        return;
    end

    % Read header
    % Read OFF identifier
    line = fgetl(fid);
    if ~strcmp(line(1:3), 'OFF');
        disp(['ERROR: file does not have the "OFF" identifier']);
        status = 1;
        return;
    end

    % Read number of vertices, faces, and edges
    line = line(4:end);
    [data, count] = sscanf(line, '%d %d %d');
    if count == 0
        line = fgetl(fid);
        [data, count] = sscanf(line, '%d %d %d');
    end
    vcount = data(1);
    fcount = data(2);
    X = zeros(vcount, 3);
    F = zeros(fcount, 3);

    % Read content
    for vnum = 1:vcount
        if feof(fid)
            disp(['ERROR: file too short']);
        end
        line = fgetl(fid);
        X(vnum, :) = sscanf(line, '%f %f %f');
    end

    for fnum = 1:fcount
        if feof(fid)
            disp(['ERROR: file too short']);
        end
        line = fgetl(fid);
        [data, count] = sscanf(line, '%d %d %d %d');
        F(fnum, :) = data(2:4) + 1;
    end

    % Close file
    fclose(fid);

    % Set up output mesh
    M = struct();
    M.vertices = X;
    M.faces = F;
