function L = norm_l1l2(x)
% L2L1 norm

    L = 0;
    for i=1:size(x,2)
        L = L + norm(x(:,i));
    end
end
