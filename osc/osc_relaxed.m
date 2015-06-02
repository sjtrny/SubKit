function [ Z ] = osc_relaxed( X, lambda_1, lambda_2, diagconstraint, mu)
%% Solves the following
%
% min || X - XZ ||_F^2 + lambda_1 || Z ||_1 + lambda_2 || Z R ||_1/2
%
% via LADMPSAP
%
% Created by Stephen Tierney
% stierney@csu.edu.au
%

if (~exist('diagconstraint','var'))
    diagconstraint = 0;
end

max_iterations = 200;

func_vals = zeros(max_iterations,1);

[~, xn, ~] = size(X);

Z = zeros(xn, xn);
Z_prev = Z;

R = (triu(ones(xn,xn-1),1) - triu(ones(xn, xn-1))) + (triu(ones(xn, xn-1),-1)-triu(ones(xn, xn-1)));
R = sparse(R);

U = zeros(xn, xn-1);
U_prev = U;

Y_2 = zeros(xn, xn-1);

if (~exist('mu','var'))
    mu = 0.1;
end

mu_max = 100;  %1;  % 0.001;
gamma_0 = 1.1;

normfX = norm(X,'fro');
rho = (norm(X)^2) * 1.1;  %1.1;
  
tol_1 = 1*10^0;
tol_2 = 1*10^-4;

covar = X'*X;
lam_2_ctrl = max(lambda_2/lambda_2, 0);

for k = 1 : max_iterations

    %Update Z
    partial = -covar  + covar*Z_prev - lam_2_ctrl*mu*(U_prev - Z_prev*R + 1/mu *Y_2)*R' ;
    
    V = Z_prev - 1/rho* partial;
    
    Z = solve_l1(V, lambda_1/rho);
    
    % Set Z diag to 0
    if (diagconstraint)
        Z(logical(eye(size(Z)))) = 0;
    end
    
    Z(Z < 0) = 0;
    
    %Update J
    partial = mu*(U_prev - Z_prev*R + 1/mu *Y_2);
    V = U_prev - 1/rho * partial;

    U = solve_l1l2(V, lambda_2/rho);
    
    Y_2 = lam_2_ctrl*(Y_2 + mu * (U - Z*R));

    % Update mu
    if (mu * max([norm(Z - Z_prev,'fro'), norm(U - U_prev,'fro')] / normfX) < tol_2)
        gamma_1 = gamma_0;
    else
        gamma_1 = 1;
    end

    mu = min(mu_max, gamma_1 * mu);
    
    % Check convergence
    func_vals(k) = .5 * norm(X - X*Z,'fro')^2 + lambda_1*norm_l1(Z) +lambda_2*norm_l1l2(Z*R);

    if ( k > 1 && norm(U - Z*R, 'fro') < tol_1 ...
            && mu * max([norm(Z - Z_prev,'fro'), norm(U - U_prev,'fro')] / normfX) < tol_2)
        break;
    end
    
    if ( k > 5 && abs(func_vals(k) - func_vals(k-1)) < 1*10^-4)
        break
    end
    
    % Update prev vals
    Z_prev = Z;
    U_prev = U;

end

end