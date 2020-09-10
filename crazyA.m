function [A,L] = crazyA(n,s)

A  = randn(n);
U1 = triu(A);
A=A+A';
b = randn(n,1);
[U V] = eig(A);

for i=1:length(b)
    if (V(i,i) < 0)
        V(i,i) = -V(i,i);
    end
end

A = U*V*inv(U) + 0.01*U1;
[P L B] = svd(A);
L(1,1) = s*L(1,1);
A = P*L*B';
A = A - 0.2*randn(n);

end