%%% TP3 : HR Methods %%%% By BOUJNOUNI Fatine 

N = 63;
f = [1/4,1/4+1/N];
a = [1,10];
delta = [0,-0.05];
phi = randn(2,1);
x = Synthesis(N,delta,f,a,phi);
Nfft = 1024;

% 3.1 - Spectral analysis by Fourier transform 

figure(1)
[pxx,w] = periodogram(x,[],N);
plot(w,10*log10(pxx),'r');
hold on;
[pxx,w] = periodogram(x,[],Nfft);
plot(w,10*log10(pxx),'b');
hold off;
legend('without zero-padding','with zero-padding');


% 3.2 High resolution methods
n = 32;
K = 2;

%%%%%%%%%%%%%%%%%%%% Function ESPRIT.m
% l = N-n+1;
% c = x(1:n);
% r = x(n:N);
% X = hankel(c,r);
% 
% Rxx = (1/l)*X*X';
% 
% [U1,Lambda,U2] = svd(Rxx);
% 
% % 3.2.1 - Esprit method 
% % 3 
% W = U1(:,1:2);
% Wd = W(1:end-1,:);
% Wu = W(2:end,:);
% 
% Phi = pinv(Wd)*Wu;
% 
% Zk = eig(Phi);
% 
% deltas = log(abs(Zk));
% freq = (1/(2*pi))*angle(Zk);

%%%%%%%%%%%%%%%%%%%%
[deltas,freq] = ESPRIT(x,n,K);

%4 
%%%%%%%%%%%%%%%%%%%% Function LeastSquares.m
% A = (0:N-1);
% B = deltas + 1j*2*pi*freq;
% C = B*A;
% 
% V = exp(transpose(C));
% 
% Alphas = pinv(V)*x;
% 
% ak = abs(Alphas);
% phi_k = angle(Alphas);
%%%%%%%%%%%%%%%%%%%%

[ak,phi_k] = LeastSquares(x,deltas,freq);


%%%%%% Comparison between synthesised values and ESPRIT+LeastSquare Values
fprintf('Frequencies synthesized:');
disp(f);
fprintf('Frequencies ESPRIT+LS:');
disp(freq');
fprintf('Deltas synthesized:');
disp(delta);
fprintf('Deltas ESPRIT+LS:');
disp(deltas);
fprintf('Amplitudes synthesized:');
disp(a);
fprintf('Amplitudes ESPRIT+LS:');
disp(ak);
fprintf('Phases synthesized:');
disp(phi);
fprintf('Phases ESPRIT+LS:');
disp(phi_k);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 3.2.2 MUSIC Method 

% W_ortho = U1(:, K+1:n);
% 
% ind_f = 1;
% ind_delta = 1;
% f_inter = 0:0.001:1; % Frequencies interval
% delta_inter = -0.1:0.001:0.1; % Deltas interval 
% for f = f_inter
%     for delta = delta_inter
%         vect = exp((delta+2*1i*pi*f)*(0:n-1)');
%         P(ind_delta, ind_f) = 1/(norm(W_ortho'*vect))^2;
%         ind_delta = ind_delta + 1;
%     end
%     ind_delta=1;
%     ind_f = ind_f+1;
% end

[f_inter, delta_inter, P] = MUSIC(x,n,K);

figure(2);
surf(f_inter, delta_inter, log10(P));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4 - Audio Signals 
[y,Fs] = audioread('ClocheA.WAV');
[y1,Fs1] = audioread('ClocheB.WAV');

F1 = fft(y);
%plot(abs(F1));

figure(3)
[pxx,w] = periodogram(y,[],N);
plot(w,10*log10(pxx),'r');
hold on;
[pxx,w] = periodogram(y,[],Nfft);
plot(w,10*log10(pxx),'b');
hold off;
legend('without zero-padding','with zero-padding');

F2 = fft(y1);
%plot(abs(F2));

figure(4)
[pxx,w] = periodogram(y1,[],N);
plot(w,10*log10(pxx),'r');
hold on;
[pxx,w] = periodogram(y1,[],Nfft);
plot(w,10*log10(pxx),'b');
hold off;
legend('without zero-padding','with zero-padding');

K = 54;
n = 512;
l = 2*n;
i = 10000;
N1 = (n+l-1)*15;

%%%%%%%%%%%%%%%%% ClocheA.WAV
x1 = y(i:(i+N1-1));

[delta_bis,f_bis] = ESPRIT(x1,n,K);

[a_bis,phi_bis] = LeastSquares(x1,delta_bis,f_bis);

x1_syn= Synthesis(N1,delta_bis,f_bis,a_bis,phi_bis);

soundsc(real(x1_syn));


%%%%%%%%%%%%%%%%% ClocheB.WAV
x2 = y1(i:(i+N1-1));

[delta_bis2,f_bis2] = ESPRIT(x2,n,K);

[a_bis2,phi_bis2] = LeastSquares(x2,delta_bis2,f_bis2);

x2_syn= Synthesis(N1,delta_bis2,f_bis2,a_bis2,phi_bis2);

soundsc(real(x2_syn));
