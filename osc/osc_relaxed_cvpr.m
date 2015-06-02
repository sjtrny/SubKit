function [ Z ] = osc_relaxed_cvpr( X, lambda_1, lambda_2, diagconstraint)
%% Solves the following
%
% min || X - XZ ||_F^2 + lambda_1 || Z ||_1 + lambda_2 || Z R ||_1/2
%
% via ADM
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

S = zeros(xn, xn); % S = Z
R = (triu(ones(xn,xn-1),1) - triu(ones(xn, xn-1))) + (triu(ones(xn, xn-1),-1)-triu(ones(xn, xn-1)));
R = sparse(R);

U = zeros(xn, xn-1);

G = zeros(xn, xn);
F = zeros(xn, xn-1); 

Z = zeros(xn, xn);

gamma_1 = 1;
gamma_2 = 1;
p = 1.1;

tol = 1*10^-3;

for k = 1 : max_iterations

    % Update Z
    V = S - (G/gamma_1);

    Z = solve_l1(V, lambda_1/gamma_1);

    % Set Z diag to 0
    if (diagconstraint)
        Z(logical(eye(size(Z)))) = 0;
    end

    % Update S
    A = X'*X + gamma_1*speye(xn,xn);
    B = gamma_2*(R*R');
    C = -(X'*X + gamma_2*U*R' + gamma_1*Z + G + F*R');

    S = lyap(A, B, C);

    % Update U
    V = S*R - (1/gamma_2)*F;

    U = solve_l1l2(V, lambda_2/gamma_2);

    % Update G, F
    G = G + gamma_1 * (Z - S);
    F = F + gamma_2 * (U - S*R);

    % Update gamma_1, gamma_2

    gamma_1 = p * gamma_1;
    gamma_2 = p * gamma_2;

    % Check convergence
    func_vals(k) = .5 * norm(X - X*Z,'fro')^2 + lambda_1*norm(Z,1) +lambda_2*norm_l1l2(Z*R);

    if k > 1
        if func_vals(k) < tol
            break
        end
    end

    if k > 100
        if func_vals(k) < tol || func_vals(k-1) == func_vals(k) ...
                || func_vals(k-1) - func_vals(k) < tol
            break
        end
    end


end

end