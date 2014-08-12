paths = ['common:', genpath('libs'), 'spatsc:'];
addpath(paths);

rng(1);

rows = 100;
n_space = 5;
cluster_size = 20;

A = rand(rows, n_space) * rand(n_space, n_space);

permute_inds = reshape(repmat(1:n_space, cluster_size, 1), 1, n_space * cluster_size );
A = A(:, permute_inds);

corruption = 0.0;

N = randn(size(A)) * corruption;

X = A + N;

X = normalize(X);

maxIterationbCluster = 200;
lambda_1 = 0.01;
lambda_2 = 0.01;
gamma_1 = 1;
gamma_2 = 1;
p = 1.1;

Z = spatsc_relaxed(X, lambda_1, lambda_2, gamma_1, gamma_2, p, maxIterationbCluster);

clusters = ncutW(abs(Z) + abs(Z'), n_space);

final_clusters = condense_clusters(clusters, 1);

imagesc(final_clusters);

rmpath(paths);