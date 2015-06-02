function [ Z ] = lrr_exact_l1( X, lambda )
%% Solves the following
%
% min || E ||_1 + lambda || Z ||_*
%   s.t. X = XZ + E
%
% Created by Stephen Tierney
% stierney@csu.edu.au
%

max_iterations = 200;

func_vals = zeros(max_iterations, 1);

Z = zeros(size(X, 2));

E = zeros(size(X));

Y = zeros(size(X));

mu = 0.1;
mu_max = 1;

normfX = norm(X,'fro');
rho = (norm(X,2)^2) * 1.2;

gamma_0 = 1.1;

tol_1 = 1*10^-2;
tol_2 = 1*10^-4;

for k = 1 : max_iterations

    Z_prev = Z;
    E_prev = E;
    
    % Solve for Z
    
    partial = mu*X'*(X*Z - (X - E - 1/mu * Y));
    
    V = Z - 1/rho * partial;
    
    [Z, s] = solve_nn(V, lambda/rho);
    
    % Solve for E
    
    V = X - X*Z - 1/mu * Y;
    
    E = solve_l1(V, 1/mu);
    
    % Update Y
    
    Y = Y + mu*(X*Z -X + E);
    
    if (mu * max(sqrt(rho) * norm(Z - Z_prev,'fro'), norm(E - E_prev))/norm(X,'fro') < tol_2)
        gamma = gamma_0;
    else
        gamma = 1;
    end
    
    mu = min(mu_max, gamma * mu);
    
    % Check convergence

    func_vals(k) = norm_l1(E) + lambda * sum(s);
    
    if ( norm(X*Z - X + E, 'fro')/normfX < tol_1 ...
            && (mu * sqrt(rho) * max([ norm(Z - Z_prev,'fro'), norm(E - E_prev, 'fro')]/normfX) < tol_2))
        break;
    end
    
end

end

