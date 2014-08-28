paths = ['common:', genpath('libs'), 'lrr:'];
addpath(paths);

rng(1);

rows = 100;
n_space = 5;
cluster_size = 20;

A = rand(rows, n_space) * rand(n_space, n_space);

permute_inds = reshape(repmat(1:n_space, cluster_size, 1), 1, n_space * cluster_size );
A = A(:, permute_inds);

corruption = 0.1;

N = randn(size(A)) * corruption;

X = A + N;

X = normalize(X);

% Z = lrr_noiseless(X, n_space);
% Z = lrr_relaxed(X, 0.01);
% Z = lrr_relaxed_lin(X, 0.01);
% Z = lrr_relaxed_lin_ext(X, 0.01);
% Z = lrr_relaxed_lin_acc(X, 0.01);
Z = lrr_exact_fro(X, 0.5);
% Z = r_lrr_fro(X, 100, n_space);

clusters = ncutW(abs(Z) + abs(Z'), n_space);

final_clusters = condense_clusters(clusters, 1);

imagesc(final_clusters);

rmpath(paths);