function [ X, s ] = solve_nn( Y, tau )
%% Solves the following
%
%   min tau * |X|_* + 1/2*|X - Y|^2
%
% Created by Stephen Tierney
% stierney@csu.edu.au
%

[U, S, V] = svd(Y, 'econ');
s = diag(S);

ind = find(s <= tau);
s(ind) = 0;

ind = find(s > tau);
s(ind) = s(ind) - tau;

S = diag(s);
    
X = U*S*(V');

end

