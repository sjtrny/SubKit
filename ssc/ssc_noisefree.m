function [ Z ] = ssc_noisefree( X )
%% Solves the following
%
% min || Z ||_1
%   s.t. X = XZ
%
% Created by Stephen Tierney
% stierney@csu.edu.au
%

sample_corr = X' * X;

sample_corr(logical(eye(size(sample_corr)))) = 0;

[~, n] = size(X);

[sorted, sorted_inds] = sort(sample_corr, 'descend');

row_inds = zeros(n^2, 1);
col_inds = zeros(n^2, 1);
count = 0;

max_vals = max(sample_corr);

for k = 1 : n

    for j = 1 : n
        
        if (abs(sorted(j,k) - max_vals(k)) < 1*10^-6)
            count = count + 1;
            
            row_inds(count) = sorted_inds(j,k);
            col_inds(count) = k;
            
        else
            break;
        end

    end

end

row_inds = row_inds(1:count);
col_inds = col_inds(1:count);

lin_inds = sub2ind(size(sample_corr), row_inds, col_inds);

Z = zeros(size(sample_corr));
Z(lin_inds) = sample_corr(lin_inds);

end

