function [ Z ] = r_lrr_fro( X, lambda, num_clusters )

A = rpca_fro(X, lambda);

Z = lrr_noiseless(A, num_clusters);

end

