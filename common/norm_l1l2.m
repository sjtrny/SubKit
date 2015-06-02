function L = norm_l1l2(x)
%% L1/L2 norm of X
%
% Created by Stephen Tierney
% stierney@csu.edu.au
%

    L = 0;
    for i=1:size(x,2)
        L = L + norm(x(:,i));
    end
end
