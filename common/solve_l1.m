function [ x ] = solve_l1( b, lambda )
%% Solves the following
%
%   arg min_{x} || x - b ||_2^2 + lambda || x ||_1
%
% Created by Stephen Tierney
% stierney@csu.edu.au
%

x = sign(b).*max(abs(b) - lambda, 0);

end

