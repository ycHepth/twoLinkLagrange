function c = Christoffel(n,D,i,j,k)
%CHRISTOFFEL just for 2-link model
%   n : freedom of coordinate
%   D : initial matrix
syms q1 q2
q = [q1 q2]';
c(i,j,k) = 1/2*(diff(D(k,j),q(i)) + diff(D(k,i),q(j)) - diff(D(i,j),q(k)));
end

