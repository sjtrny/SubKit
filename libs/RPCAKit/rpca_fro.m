function [ A ] = rpca_fro( X, lambda )
%% Solves the following
%
% min || A ||_* + lambda/2 * || N ||_F^2
%   s.t. X = A + N
%
% Created by Stephen Tierney
% stierney@csu.edu.au
%

max_iterations = 100;

func_vals = zeros(max_iterations, 1);

A = zeros(size(X));

N = zeros(size(X));

Y = zeros(size(X));

mu = 0.5;
mu_max = 100;

gamma_0 = 1.1;

normfX = norm(X,'fro');

tol_1 = 1*10^-2;
tol_2 = 1*10^-4;

for k = 1 : max_iterations

    A_prev = A;
    N_prev = N;
    
    % Solve for A
    
    V = X - N + 1/mu * Y;
    
    [A, s] = solve_nn(V, 1/mu);
    
    % Solve for N
    
    V = X - A + 1/mu * Y;
    
    N = solve_l2(V, lambda/mu);
    
    % Update Y
    
    Y = Y + mu*(X - A - N);
    
    if (mu * max(norm(A - A_prev,'fro'), norm(N - N_prev))/norm(X,'fro') < 1*10^-3)
        gamma = gamma_0;
    else
        gamma = 1;
    end
    
    mu = min(mu_max, gamma * mu);
    
    % Check convergence
    
    func_vals(k) = sum(s) + lambda*0.5*norm(N, 'fro')^2;

    if ( norm(X - A - N, 'fro') < tol_1 ...
            && (mu * max([ norm(A - A_prev,'fro'), norm(N - N_prev, 'fro')]) / normfX < tol_2))
        break;
    end
    
    
end

end