function [ coeff ] = gen_depcoeff( m, v, dim_space, cluster_size )
%% Generate dependant coefficients
%
% Created by Stephen Tierney
% stierney@csu.edu.au
%
    v_nbor = v / 2;

    lower_sigma = v_nbor*(triu(ones(cluster_size, cluster_size),-1) - triu(ones(cluster_size, cluster_size)));

    sigma = tril(lower_sigma)' + tril(lower_sigma);
    sigma(logical(eye(size(sigma)))) = v*ones(cluster_size, 1);

    coeff = zeros(dim_space, cluster_size);
    for k = 1 : dim_space
        coeff(k, :) = mvnrnd( m + zeros(cluster_size, 1), sigma , 1);
    end


end

