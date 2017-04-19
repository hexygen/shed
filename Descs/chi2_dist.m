function [ d ] = chi2_dist( hist1, hist2 )
% CHI2_DIST Computes chi-squared distance between the histograms
%
% Written by Yanir Kleiman, 2014

s = hist1 + hist2;

diff = (hist1 - hist2).^2 ./ s;
diff(s == 0) = 0;
d = sum(diff);


end

