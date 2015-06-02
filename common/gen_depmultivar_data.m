function [ X ] = gen_depmultivar_data( dim_data, dim_space, cluster_size, n_space, m, v )
%% Generate Dependant Multivariate Data
%
% Created by Stephen Tierney
% stierney@csu.edu.au
%

    X = zeros(dim_data, cluster_size*n_space);

    [T, ~] = qr(rand(dim_data,dim_data));
    if det(T) < 1
        T(:,1) = -T(:,1);
    end
    
    basis = orth(randn(dim_data, dim_space));

    X(:, 1:cluster_size) = basis * gen_depcoeff(m, v, dim_space, cluster_size);

    for j = 1 : n_space-1
        
        basis = T * basis;
                
        X(:, cluster_size*j+1: cluster_size*(j+1)) = basis * gen_depcoeff(m, v, dim_space, cluster_size);
    
    end


end

