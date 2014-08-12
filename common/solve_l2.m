function [ x ] = solve_l2( w, lambda )
%MYSOLVE_L2 Summary of this function goes here
%   Detailed explanation goes here

% min lambda \lambda |x|_2 + 1/2 * |x-w|_2^2

% x = w ./ (2*lambda + 1);

x = w ./ (lambda + 1);

end

