paths = ['common:', genpath('libs'), 'osc:'];
addpath(paths);

rng(1);

rows = 100;
n_space = 5;
cluster_size = 50;

A = rand(rows, n_space) * rand(n_space, n_space);

permute_inds = reshape(repmat(1:n_space, cluster_size, 1), 1, n_space * cluster_size );
A = A(:, permute_inds);

corruption = 0;

N = randn(size(A)) * corruption;

X = A + N;

X = normalize(X);

% lambda_1 = 0.00000001;
% lambda_2 = 0.00000001;
% 
% Z = osc_relaxed_cvpr(X, lambda_1, lambda_2);

lambda_1 = 0.099;
lambda_2 = 0.001;
Z = osc_relaxed(X, lambda_1, lambda_2);

% lambda_1 = 0.99;
% lambda_2 = 0.001;
% Z = osc_exact(X, lambda_1, lambda_2);

clusters = ncutW(abs(Z) + abs(Z'), n_space);

final_clusters = condense_clusters(clusters, 1);

imagesc(final_clusters);

rmpath(paths);