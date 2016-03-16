function [ C ] = contr( A, B )
%CONTR contraction of two tensors
%   Detailed explanation goes here

    %C = trace(A*B');
    
    %for i=1:size(A,1)
    %   C(i) = sum(A(i,:)*B(i,:)); % should be B(:,i) ?
    %end

    C = zeros(size(A,3),1);
    for i=1:size(C,1)
       C(i) = sum(sum(A(:,:,i).*B(:,:)'));
    end

    
end

