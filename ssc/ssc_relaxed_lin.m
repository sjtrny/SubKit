function [ Z ] = ssc_relaxed_lin( X, lambda )
%% Solves the following
%
% min || X - XZ ||_F^2 + lambda || Z ||_1
%
% by gradient descent
%
% Created by Stephen Tierney
% stierney@csu.edu.au
%

max_iterations = 100;

func_vals = zeros(max_iterations, 1);
previous_func_val = Inf;

n = size(X, 2);

Z = zeros(n);

rho = norm(X,2)^2;

for k = 1 : max_iterations

    % Solve for Z
    
    partial = -X'*X + X'*X*Z;
    
    V = Z - 1/rho * partial;
    
    Z = solve_l1(V, lambda/rho);
    
    Z(logical(eye(size(Z)))) = 0;

    % Check convergence

    func_vals(k, 1) = 0.5*norm(X - X*Z, 'fro')^2 + lambda*norm_l1(Z);
    
    if ( abs(func_vals(k) - previous_func_val) <= 1*10^-6 )
        break;
    else
        previous_func_val = func_vals(k);
    end
   
end

end

