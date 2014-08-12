function [ Z ] = r_lrr_l1l2( X, lambda, num_clusters )

A = rpca_l1l2(X, lambda);

Z = lrr_noiseless(A, num_clusters);

end

