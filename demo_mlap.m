paths = ['common:', genpath('libs'), 'mlap:'];
addpath(paths);

rng(1)

n_spectra = 5;

n_cluster = 5;
n_space = n_cluster;
cluster_size = 20;
m = 0.1;
v = 0.001;

truth = reshape(repmat(1:n_cluster,cluster_size,1),1,n_cluster*cluster_size)';

B1 = gen_depmultivar_data(200, 4, cluster_size, n_space, m, v);
B2 = gen_depmultivar_data(50, 6, cluster_size, n_space, m, v);
B3 = gen_depmultivar_data(100, 5, cluster_size, n_space, m, v);

noise_mag = 0;

noise = randn(size(B1));
w = noise * noise_mag;
X1 = B1 + w;

noise = randn(size(B2));
w = noise * noise_mag;
X2 = B2 + w;

noise = randn(size(B3));
w = noise * noise_mag;
X3 = B3 + w;

X_1normed = normalize(X1);
X_2normed = normalize(X2);
X_3normed = normalize(X3);

Xs = {X_1normed, X_2normed, X_3normed};

alpha = 1;
lambda = 0.01;
[Zs] = mlap(Xs, alpha, lambda);

Z_final = sqrt(sum(Zs.^2, 3));

clusters = ncutW((abs(Z_final)+abs(Z_final')), n_space);

final_clusters = condense_clusters(clusters, 1);
