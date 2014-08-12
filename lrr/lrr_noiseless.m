function [ Z ] = lrr_noiseless( A, num_clusters )
%LRR_NOISELESS Summary of this function goes here
%   Detailed explanation goes here

[~, ~, V] = svd(A);

Z = V(:,1:num_clusters) * V(:,1:num_clusters)';

end

