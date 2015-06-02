function [ Z ] = lsr_relaxed( X, lambda, diag )
%% Solves the following
%
% min 1/2 || X - XZ ||_F^2 + lambda/2 || Z ||_F^2
%   s.t. diag(Z) = 0 (when diag = 1)
%
% Created by Stephen Tierney
% stierney@csu.edu.au
%

if (diag)
    D = (X' * X + lambda*eye(size(X,2)));
    Z = spdiags(diag(D), 0, size(X,2), size(X,2)) \ D;
    Z(logical(eye(size(Z)))) = 0;
else
    Z = (X' * X + lambda*eye(size(X,2))) \ X'*X;
end

end

