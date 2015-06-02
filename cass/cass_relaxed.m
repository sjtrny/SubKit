function [ Z ] = cass_relaxed( X, lambda )
%% Solves the following
%
% min 1/2 || X - XZ ||_F^2 + \sum_i^N lambda || X diag(Z(:,i) ||_*
%
% via ADMAP. Since each col of Z is independant we handle them one at a time.
%
% Created by Stephen Tierney
% stierney@csu.edu.au
%

max_iterations = 200;

N = size(X, 2);

Z = zeros(N);

mu_max = 1*10^1;

gamma_0 = 1.1;

tol_1 = 1*10^-2;
tol_2 = 1*10^-4;

for j = 1 : N
    
    func_vals = zeros(max_iterations, 1);
    
    mu = 0.1;
    
    z = zeros(N-1, 1);
    J = zeros(size(X,1), N-1);
    Y = zeros(size(X,1), N-1);
    
    sub_X = X(:, 1:end ~= j);
    
    XtX = sub_X' * sub_X;
    diagXtX = diag(diag(XtX));
    
    for k = 1 : max_iterations
        
        z_prev = z;
        J_prev = J;
        
        % Update J
        [J, s] = solve_nn(sub_X * diag(z) - 1/mu * Y, lambda/mu);
        
        % Update z
        z = (XtX + mu * diagXtX) \ (sub_X' * X(:,j) + diag(sub_X'*(Y + mu * J)) );
        
        % Update Y
        Y = Y + mu * (J - sub_X * diag(z));
        
        if max(max(abs(z - z_prev)), max(max(abs(J - J_prev))) ) > tol_2
            gamma = gamma_0 ;
        else
            gamma = 1 ;
        end
        
        mu = min(mu_max, gamma * mu);
    
        % Check convergence

        func_vals(k) = 0.5*norm(X(:,j) - sub_X * z, 'fro')^2 + lambda*sum(s);

        if ( max(max(abs(J - sub_X * diag(z)))) < tol_1...
            && max(max(abs(J - J_prev))) < tol_1...
            && max(max(abs(z - z_prev))) < tol_1 )
            break;
        end
        
    end
    
    Z(1:end ~= j,j) = z;
end


end

