function [delta,f] = ESPRIT(x,n,K)

N = size(x,1);
l = N-n+1;
c = x(1:n);
r = x(n:N);
X = hankel(c,r);

Rxx = (1/l)*X*X';

[U1,Lambda,U2] = svd(Rxx);

W = U1(:,1:K);
Wd = W(1:end-1,:);
Wu = W(2:end,:);

Phi = pinv(Wd)*Wu;

Zk = eig(Phi);

delta = log(abs(Zk));
f = (1/(2*pi))*angle(Zk);
end