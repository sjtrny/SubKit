function [ E ] = solve_l1l2( W, lambda )
%% Solves the following
%
%   arg min_{x} || X - W ||_2^2 + lambda || X ||_1/2
%
% Created by Stephen Tierney
% stierney@csu.edu.au
%

[m, n] = size(W);

E = W;

for i = 1 : n
    
    norm_col = norm(W(:,i));
    
    if (norm_col > lambda)
        E(:,i) = (norm_col - lambda) * W(:,i) / norm_col;
    else
        E(:,i) = zeros(m, 1);
    end
    
end

end
