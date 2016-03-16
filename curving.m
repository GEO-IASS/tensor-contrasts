function [ C ] = curving( dphi2, dphi3, D )
%CURVING Calculates the curving index
%   Detailed explanation goes here

[E,L] = eig(D);
L = diag(L);
[~,I] = sort(L,1,'descend');
e1 = E(:,I(1));

C = sqrt( abs(contr(tprod(dphi2,dphi2) + tprod(dphi3,dphi3), tprod(e1,e1)) ));

end

