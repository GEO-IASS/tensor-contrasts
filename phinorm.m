function [ phi ] = phinorm( p, D )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[E,L] = eig(D);
L = diag(L);
[~,I] = sort(L,1,'descend');

e1 = E(:,I(1));
e2 = E(:,I(2));
e3 = E(:,I(3));

phi = [];
switch p
    case 1
        phi = (tprod(e2, e3) + tprod(e3, e2))/sqrt(2);
        
    case 2
        phi = (tprod(e3, e1) + tprod(e1, e3))/sqrt(2);
        
    case 3
        phi = (tprod(e1, e2) + tprod(e2, e1))/sqrt(2);
end


end

