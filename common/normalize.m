function x = normalize(y)
%% Normalizes each column such that the norm of each col is 1
%
% Created by Stephen Tierney
% stierney@csu.edu.au
%

    for i=1:size(y,2)
        x(:,i) = y(:,i)/sqrt(y(:,i)'*y(:,i));
    end
end