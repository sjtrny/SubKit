function [ Z ] = spatsc_relaxed( X, lambda_1, lambda_2, diagconstraint)

if (~exist('diagconstraint','var'))
    diagconstraint = 0;
end

max_iterations = 30;

func_vals = zeros(max_iterations,1);

[~, xn, ~] = size(X);

Z = zeros(xn, xn);
Z_prev = Z;

S = zeros(xn, xn);
S_prev = S;

R = (triu(ones(xn,xn-1),1) - triu(ones(xn, xn-1))) + (triu(ones(xn, xn-1),-1)-triu(ones(xn, xn-1)));
R = sparse(R);

U = zeros(xn, xn-1);
U_prev = U;

Y_1 = zeros(xn, xn);
Y_2 = zeros(xn, xn-1);

normfX = norm(X,'fro');

mu = 0.0001;
mu_max = 10;
gamma_0 = 1.1;

lambda_1 = lambda_1*mu;

tol_1 = 1*10^-2;
tol_2 = 1*10^-2;

for k = 1 : max_iterations

    % Update Z
    V = S_prev - (Y_1 / mu);

    Z = solve_l1(V, lambda_1/mu);

    % Set Z diag to 0
    if (diagconstraint)
        Z(logical(eye(size(Z)))) = 0;
    end

    % Update S
    A = X'*X + mu*speye(xn,xn);
    B = mu*(R*R');
    C = -(X'*X + mu*U_prev*R' + mu*Z_prev + Y_1 + Y_2*R');

    S = lyap(A, B, C);

    % Update U
    V = S_prev*R - (Y_2 / mu);

    U = solve_l1(V, lambda_2/mu);

    % Update Y_1 and Y_2
    Y_1 = Y_1 + mu * (Z - S);
    Y_2 = Y_2 + mu * (U - S*R);

    % Update mu
    if (mu * max([norm(Z - Z_prev,'fro'), norm(S - S_prev), norm(U - U_prev,'fro'), norm(S*R - S_prev*R,'fro')] / normfX) < tol_2)
        gamma_1 = gamma_0;
    else
        gamma_1 = 1;
    end
    
    mu = min(mu_max, gamma_1 * mu);
    
    % Check convergence
    func_vals(k) = .5 * norm(X - X*Z,'fro')^2 + lambda_1*norm_l1(Z) +lambda_2*norm_l1(Z*R);

    if ( norm(Z - S, 'fro')/ normfX < tol_1 && norm(U - S*R, 'fro')/ normfX < tol_1 ...
            && mu * max([norm(Z - Z_prev,'fro'), norm(S - S_prev), norm(U - U_prev,'fro'), norm(S*R - S_prev*R,'fro')] / normfX) < tol_2)
        break;
    end
    
    if ( k > 5 && abs(func_vals(k) - func_vals(k-1)) < 1*10^-3)
        break
    end
    
    % Update prev vals
    Z_prev = Z;
    S_prev = S;
    U_prev = U;

end

end