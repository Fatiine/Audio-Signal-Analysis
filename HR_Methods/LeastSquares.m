function [a,phi] = LeastSquares(x,delta,f)
N = size(x,1);

A = (0:N-1);
B = delta + 1j*2*pi*f;
C = B*A;

V = exp(transpose(C));

Alphas = pinv(V)*x;

a = abs(Alphas);
phi = angle(Alphas);

end