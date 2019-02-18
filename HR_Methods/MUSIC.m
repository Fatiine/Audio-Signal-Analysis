function [f_vec, delta_vec, P] = MUSIC(x,n,K)
N = size(x,1);
l = N-n+1;
c = x(1:n);
r = x(n:N);
X = hankel(c,r);
Rxx = (1/l)*X*X';

[U1,Lambda,U2] = svd(Rxx);

W_ortho = U1(:, K+1:n);

ind_f = 1;
ind_delta = 1;
f_vec = 0:0.001:1; % Frequencies interval
delta_vec = -0.1:0.001:0.1; % Deltas interval 
for f = f_vec
    for delta = delta_vec
        vec = exp((delta+2*1i*pi*f)*(0:n-1)');
        P(ind_delta, ind_f) = 1/(norm(W_ortho'*vec))^2;
        ind_delta = ind_delta + 1;
    end
    ind_delta=1;
    ind_f = ind_f+1;
end

end