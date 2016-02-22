function [ x ] = solve_l2( w, lambda )
%% Solves the following
%
%   arg min_{x} 1/2 || x - b ||_2^2 + lambda/2 || x ||_2^2
%
% Created by Stephen Tierney
% stierney@csu.edu.au
%

x = w ./ (lambda + 1);

end

