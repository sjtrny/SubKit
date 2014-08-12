function [ val ] = norm_l1( X )
%NORM_L1 Summary of this function goes here
%   Detailed explanation goes here

val = sum(sum(abs(X)));

end

