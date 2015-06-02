function [ Z ] = lrr_noiseless( A, num_clusters )
%% Solves the following
%
% min || Z ||_*
%   s.t. X = XZ
%
% Created by Stephen Tierney
% stierney@csu.edu.au
%

[~, ~, V] = svd(A);

Z = V(:,1:num_clusters) * V(:,1:num_clusters)';

end

