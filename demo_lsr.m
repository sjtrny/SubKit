paths = ['common:', genpath('libs'), 'lsr:'];
addpath(paths);

rng(1);

dim_data = 100;
dim_space = 4;
n_space = 5;
cluster_size = 20;
m = 0.1;
v = 0.001;

A = gen_depmultivar_data(dim_data, dim_space, cluster_size, n_space, m, v);
A = normalize(A);

corruption = 0;

N = randn(size(A)) * corruption;

X = A + N;

X = normalize(X);

Z = lsr_relaxed(X, 0.01);

clusters = ncutW(abs(Z) + abs(Z'), n_space);

final_clusters = condense_clusters(clusters, 1);

figure, imagesc(final_clusters);

rmpath(paths);