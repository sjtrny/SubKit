function [ Z ] = r_lrr_fro( X, lambda, num_clusters )
%% Solves the following
%
% min || Z ||_*
%   s.t. A = AZ, X = A + N
%
%   Where N is assumed to be Gaussian
%
%   This problem is decomposable into RPCA then Noiseless LRR
%   lambda is therefore the paramter for RPCA
%
% Created by Stephen Tierney
% stierney@csu.edu.au
%

A = rpca_fro(X, lambda);

Z = lrr_noiseless(A, num_clusters);

end

