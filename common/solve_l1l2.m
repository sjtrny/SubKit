function [ E ] = solve_l1l2( W, lambda )
%MYSOLVE_L1L2 Summary of this function goes here
%   Detailed explanation goes here

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
