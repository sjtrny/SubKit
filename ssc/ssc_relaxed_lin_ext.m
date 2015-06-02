function [ Z ] = ssc_relaxed_lin_ext( X, lambda )
%% Solves the following
%
% min || X - XZ ||_F^2 + lambda || Z ||_1
%
% by extended gradient descent
%
% Created by Stephen Tierney
% stierney@csu.edu.au
%

max_iterations = 100;

func_vals = zeros(max_iterations, 1);
previous_func_val = Inf;

n = size(X, 2);

Z = zeros(n);

rho = 1;

gamma = 1.1;

for k = 1 : max_iterations

    Z_prev = Z;
    
    searching = true;
    while( searching )
        partial = -X'*X + X'*X*Z_prev;

        V = Z_prev - 1/rho * partial;

        Z = solve_l1(V, lambda/rho);
        
        Z(logical(eye(size(Z)))) = 0;
    
        func_vals(k, 1) =  lambda*norm_l1(Z) + 1/2 *norm(X - X*Z, 'fro')^2;
        
        approx = 0.5*norm(X - X*Z_prev, 'fro')^2 + sum(sum((Z - Z_prev).*partial)) + 0.5*rho*norm(Z - Z_prev, 'fro')^2 +  lambda*norm_l1(Z);
  
        if ( func_vals(k, 1) > approx )
            rho = gamma * rho;
        else
            searching = false;
        end
    end
    
    if ( abs(func_vals(k) - previous_func_val) <= 1*10^-6 )
        break;
    else
        previous_func_val = func_vals(k);
    end

end

end

