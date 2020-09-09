function Cq = TwoLinkDirectC(n,D,q,dq)
%TWOLINKDIRECTC Summary of this function goes here
%   C matrix of lagrange equation


for k = 1:n
    for j = 1:n
        for i = 1:n
            Ct(i) = 1/2*((diff(D(k,j),q(i)) + diff(D(k,i),q(j)) - diff(D(i,j),q(k))))*dq(i);
        end
        C(k,j) = sum(Ct);
    end
end

Cq = simplify(C);
end
