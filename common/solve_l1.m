function [ x ] = solve_l1( b, lambda )
% Function to solve soft thresholding problem
%
% arg min_{x} ||x - b||_{2}^{2} + lambda*||x||_{1}
%

x = zeros(size(b));

k = find(b > lambda);
x(k) = b(k) - lambda; 

k = find(abs(b) <= lambda);
x(k) = 0; 

k = find(b < - lambda);
x(k) = b(k) + lambda;

end

