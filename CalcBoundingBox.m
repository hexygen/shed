function [ BB ] = CalcBoundingBox( points )
%%% CALCBOUNDINGBOX Calculates a PCA-based bounding box of the part
% and returns its position, rotation, and scale.
%
% Input: a set of points that represents the shape or segment.
%
% Output:
%   - BB.rot is the rotation of the bounding box
%   - BB.scale is an *estimation* of the scale of the bounding box (this
%   is not the exact scale!)
%   - BB.volume is an *estimation* of the volume of the shape (this is not
%   the mathematical volume of the shape!)
%
%
%%% If you use this code, please cite the following paper:
%  
%  SHED: Shape Edit Distance for Fine-grained Shape Similarity 
%  Yanir Kleiman, Oliver van Kaick, Olga Sorkine-Hornung, Daniel Cohen-Or 
%  SIGGRAPH ASIA 2015
%
%%% Copyright (c) 2015 Yanir Kleiman <yanirk@gmail.com>


n = size(points, 1);

BB.center = mean(points);

% Normalize points:
pn = points - repmat(BB.center, n, 1);

% Finding the best rotation of the shape according to PCA (svd):
% (normalizing the scale by dividing the svd by the number of samples)
[U, S, V] = svd(pn'*pn / n);
BB.rot = U;

% Returning the standard deviation on each axis - this can be changed to 
% max / some percentile / some other approximation of the scale:
BB.scale = sqrt(diag(S));


% This is an estimation of "shape volume" which I found to reflect more or
% less human intuition about the volume of a shape. A shape which is tall
% and narrow tends to look like it has less volume than it truly has, thus
% using the real volume (prod(BB.scale)) sometimes give unexpected results.
% The following estimation gives more "weight" to longer sides of the
% shape, so if the longest side changes it affects the "volume" more than
% if the shortest side changes.
alpha = 0.75;

BB.volume = sum(BB.scale .^ alpha);

end

