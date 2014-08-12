function [ Z, funVal, iteration ] = osc_relaxed( X, lambda_1, lambda_2, gamma_1, gamma_2, p, maxIterations, diagconstraint)

if (~exist('diagconstraint','var'))
    diagconstraint = 0;
end

funVal = zeros(maxIterations,1);

[~, xn, ~] = size(X);

S = zeros(xn, xn); % S = Z
R = (triu(ones(xn,xn-1),1) - triu(ones(xn, xn-1))) + (triu(ones(xn, xn-1),-1)-triu(ones(xn, xn-1)));
R = sparse(R);

U = zeros(xn, xn-1);

G = zeros(xn, xn);
F = zeros(xn, xn-1); 

Z = zeros(xn, xn);

for iteration=1:maxIterations

    %% Step 1
    V = S - (G/gamma_1);
    [vm, vn, ~] = size(V);
    rolled_v = reshape(V,vm*vn,1);

    rolled_z = solve_l1(rolled_v, lambda_1/gamma_1);

    Z = reshape(rolled_z, vm, vn);

    % Set Z diag to 0
    if (diagconstraint)
        Z(logical(eye(size(Z)))) = 0;
    end

    %% Step 2
    A = X'*X + gamma_1*speye(xn,xn);
    B = gamma_2*(R*R');
    C = -(X'*X + gamma_2*U*R' + gamma_1*Z + G + F*R');

    S = lyap(A, B, C);

    %% Step 3
    V = S*R - (1/gamma_2)*F;

    U = solve_l1l2(V, lambda_2/gamma_2);

    %% Step 4

    G = G + gamma_1 * (Z - S);

    %% Step 5

    F = F + gamma_2 * (U - S*R);

    %% Step 6

    gamma_1 = p * gamma_1;
    gamma_2 = p * gamma_2;

    %% Calculate function values
    funVal(iteration) = .5 * norm(X - X*Z,'fro')^2 + lambda_1*norm(Z,1) +lambda_2*norm_l1l2(Z*R);

    if iteration > 1
        if funVal(iteration) < 1*10^-3
            break
        end
    end

    if iteration > 100
        if funVal(iteration) < 1*10^-3 || funVal(iteration-1) == funVal(iteration) ...
                || funVal(iteration-1) - funVal(iteration) < 1*10^-3
            break
        end
    end


end

end