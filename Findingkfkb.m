
n = 100;
[A,S] = crazyA(n,10^4);
%A = randn(n);
% [U S V] = svd(A);
% S = diag(linspace(1,10^6,n));
% A = U*S*V';
NORM_A = norm(A);
COND = cond(A);

for i = 1:10000
    b = randn(n,1);
    xt = A\b;
    kf(i) = NORM_A*norm(xt)/norm(b);
    kb(i) = COND/kf(i);
end

plot(1:10000,sort(kf),'r');
hold on;
plot(1:10000,sort(kb),'g');