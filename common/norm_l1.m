function [ val ] = norm_l1( X )
%% L1 norm of X
%
% Created by Stephen Tierney
% stierney@csu.edu.au
%

val = sum(abs(X(:)));

end

