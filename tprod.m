function [ C ] = tprod( A, B )
%TPROD tensor product
%   Detailed explanation goes here
    C = kron(A,B);
    C = reshape(C,3,3);
end

