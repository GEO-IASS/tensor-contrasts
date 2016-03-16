function [ D ] = dispersion( dphi2, dphi3, D )
%DISPERSION Calculates the dispersion index
%   Detailed explanation goes here

[E,L] = eig(D);
L = diag(L);
[~,I] = sort(L,1,'descend');
e2 = E(:,I(2));
e3 = E(:,I(3));

D = sqrt( abs(contr(tprod(dphi2,dphi2) + tprod(dphi3,dphi3), (tprod(e2,e2) + tprod(e3,e3))) ));

end

