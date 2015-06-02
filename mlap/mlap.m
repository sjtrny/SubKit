function [ Zs ] = mlap( Xs, lambda, tau )
%% Solves the following
%
% min alpha || A ||_1/2 + \sum_i^K (lambda || E_i ||_1/2 +  || Z_i ||_*)
%   s.t. X_i = X_i Z_i + E_i
%        A = [Z_1; Z_2; ... ; Z_K]
%
% Created by Stephen Tierney
% stierney@csu.edu.au
%

if (length(lambda) == 1)
    lambda = repmat(lambda, size(Xs, 2), 1);
end

max_iterations = 200;

func_vals = zeros(max_iterations, 1);

num_obs = length(Xs);
N = size(Xs{1}, 2);


Zs = zeros(N, N, num_obs);

Es = cell(size(Xs));
Ys = cell(size(Xs));

for k = 1 : num_obs
    Es{k} = zeros(size(Xs{k}));
    Ys{k} = zeros(size(Xs{k}));
end

% A is tall, each column corresponds to a Z
A = zeros(N^2, num_obs);
Ws = zeros(N, N, num_obs);

mu = 0.1;
mu_max = 1;

rho_list = zeros(num_obs, 1);
for k = 1 : num_obs
    rho_list(k) = (norm(Xs{k},2)^2) * 1.1;
end
rho = max(rho_list);

gamma_0 = 1.1;

tol_1 = 1*10^-2;
tol_2 = 1*10^-4;

for t = 1 : max_iterations

    Zs_prev = Zs;
    Es_prev = Es;
    A_prev = A;

    sub_func_val = 0;
    
    for k = 1 : num_obs
        % Solve for Z

        partial = mu*Xs{k}'*(Xs{k}*Zs_prev(:,:,k) - (Xs{k} - Es_prev{k} - 1/mu * Ys{k}))...
            + mu*(Zs_prev(:,:,k) - reshape(A_prev(:,k), size(Zs(:,:,k))) + 1/mu * Ws(:,:,k));

        V_Z = Zs_prev(:,:,k) - 1/rho * partial;

        [Zs(:,:,k), s_z] = solve_nn(V_Z, lambda(k)/rho);
        
        % Solve for E

        V_E = Xs{k} - Xs{k}*Zs_prev(:,:,k) - 1/mu * Ys{k};
        Es{k} = solve_l2(V_E, 1/mu);

        % Update Y

        Ys{k} = Ys{k} + mu*(Xs{k}*Zs(:,:,k) - Xs{k} + Es{k});
        sub_func_val = sub_func_val + lambda(k)*sum(s_z) + 0.5*norm(Es{k}, 'fro')^2;
        
    end

    % Solve for A
    
    V_A = reshape(Zs_prev + 1/mu * Ws, N^2, num_obs);
    
    [A, s_a] = solve_nn(V_A, tau/mu);

    % Update W
    
    Ws = Ws + reshape(mu*(reshape(Zs_prev, N^2, num_obs) - A), [N, N, num_obs]);
    
    % Update function value
    
    func_vals(t) = sub_func_val + tau*sum(s_a);

    % Check convergence
    
    fit_error = zeros(num_obs, 1);
    az_error = zeros(num_obs, 1);
    Z_shp = reshape(Zs, N^2, num_obs);
    z_diff = zeros(num_obs, 1);
    e_diff = zeros(num_obs, 1);
    for k = 1: num_obs
        fit_error(k) =  norm(Xs{k}*Zs(:,:,k) - Xs{k} + Es{k}, 'fro')/ norm(Xs{k}, 'fro');
        az_error(k) = norm(A(:,k) - Z_shp(:,k), 2);
        z_diff(k) = norm(Zs(:,:,k) - Zs_prev(:,:,k), 'fro')/ norm(Xs{k}, 'fro');
        e_diff(k) = norm(Es{k} - Es_prev{k}, 'fro')/ norm(Xs{k}, 'fro');
    end
    
    if ( ( max(fit_error) < tol_1 && max(az_error) < tol_1) ...
            && (mu * sqrt(rho)* max([ max(z_diff), max(e_diff), norm(A-A_prev,'fro')]) < tol_2))
        break;
    end
    
    % Update mu
    if (mu *sqrt(rho)* max([ max(z_diff), max(e_diff), norm(A-A_prev,'fro')]) < tol_2)
        gamma = gamma_0;
    else
        gamma = 1;
    end
     
    mu = min(mu_max, gamma * mu);  
    

end

end

