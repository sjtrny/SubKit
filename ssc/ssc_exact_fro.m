function [ Z ] = ssc_exact_fro( X, lambda )

max_iterations = 200;

func_vals = zeros(max_iterations, 1);
previous_func_val = Inf;

Z = zeros(size(X, 2));

E = zeros(size(X));

Y = zeros(size(X));

mu = 0.1;
mu_max = 1;
rho = (norm(X,2)^2) * 1.2;

gamma_0 = 1.1;

for k = 1 : max_iterations

    Z_prev = Z;
    E_prev = E;
    
    % Solve for Z
    
    partial = mu*X'*(X*Z - (X - E - 1/mu * Y));
    
    V = Z - 1/rho * partial;
    
    Z = solve_l1(V, lambda/rho);
    
    Z(logical(eye(size(Z)))) = 0;
    
    % Solve for E
    
    V = X - X*Z - 1/mu * Y;
    
    E = solve_l2(V, 1/mu);
    
    % Update Y
    
    Y = Y + mu*(X*Z -X + E);
    
    if (mu * max(sqrt(rho) * norm(Z - Z_prev,'fro'), norm(E - E_prev))/norm(X,'fro') < 1*10^-3)
        gamma = gamma_0;
    else
        gamma = 1;
    end
    
    mu = min(mu_max, gamma * mu);
    
    % Check convergence
    
    func_vals(k) = 0.5*norm(E, 'fro')^2 + lambda*norm_l1(Z);
    
    if ( abs(func_vals(k) - previous_func_val) <= 1*10^-6 )
        break;
    else
        previous_func_val = func_vals(k);
    end
    
end

end

