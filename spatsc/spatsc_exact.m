function [ Z ] = spatsc_exact( X, lambda_1, lambda_2, diagconstraint )
%% Solves the followingr
%
% min || E ||_F^2 + lambda_1 || Z ||_1 + lambda_2 || Z R ||_1
%   s.t. X = XZ + E
%
% Created by Stephen Tierney
% stierney@csu.edu.au
%

if (~exist('diagconstraint','var'))
    diagconstraint = 0;
end

max_iterations = 200;

func_vals = zeros(max_iterations,1);

[xm, xn, ~] = size(X);

Z = zeros(xn, xn);
Z_prev = Z;

E = zeros(xm, xn);
E_prev = E;

R = (triu(ones(xn,xn-1),1) - triu(ones(xn, xn-1))) + (triu(ones(xn, xn-1),-1)-triu(ones(xn, xn-1)));
R = sparse(R);

J = zeros(xn, xn-1);
J_prev = J;

Y_1 = zeros(xm, xn);
Y_2 = zeros(xn, xn-1);

if (~exist('mu','var'))
    mu = 0.1;
end

mu_max = 10;  %1;  % 0.001;
gamma_0 = 1.1;

normfX = norm(X,'fro');
rho = (norm(X)^2) * 1.1;  %1.1;
  
tol_1 = 1*10^-2;
tol_2 = 1*10^-4;

for k = 1 : max_iterations

    % Update Z
    partial =  mu*(X' * (X * Z_prev - (X - E_prev - 1/mu * Y_1)) + (Z_prev*R - (J_prev + 1/mu * Y_2)) * R');
    
    V = Z_prev  -  1/rho* partial;
    
    Z = solve_l1(V, lambda_1/(rho));

    % Set Z diag to 0
    if (diagconstraint)
        Z(logical(eye(size(Z)))) = 0;
    end
    
    % Update E
    
    V = -X*Z_prev + X - 1/mu * Y_1;
    
    E = solve_l2(V, 1/mu);
    
    % Update J    
    partial = mu*(J_prev - Z_prev*R + 1/mu *Y_2);
    V = J_prev - 1/rho * partial;

    J = solve_l1(V, lambda_2/rho);

    % Update Y_1 and Y_2
    Y_1 = Y_1 + mu * (X*Z - X + E);
    Y_2 = Y_2 + mu * (J - Z*R);
    
    % Update mu
    if (mu * (max([sqrt(rho)*norm(Z - Z_prev,'fro'), norm(E - E_prev, 'fro'), norm(J-J_prev,'fro')] / normfX)) < tol_2)
        gamma = gamma_0;
    else
        gamma = 1;
    end
     
    mu = min(mu_max, gamma * mu);    

    func_vals(k) = .5 * norm(E,'fro')^2 + lambda_1*norm_l1(Z) +lambda_2*norm_l1(Z*R);
    
    % Check convergence
    if ( (norm(X*Z - X + E, 'fro')/normfX < tol_1 && norm(J - Z*R, 'fro')/normfX < tol_1) ...
            && (mu * max([ sqrt(rho)*norm(Z - Z_prev,'fro'), norm(E - E_prev, 'fro'), norm(J-J_prev,'fro')]) / normfX < tol_2))
        break;
    end
     
    % Update prev vals
    Z_prev = Z; 
    E_prev = E;
    J_prev = J; 
    
end

end
