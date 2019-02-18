%Script for multi-pitch detection (reference: Klapuri)
% but with pitch detection with spectral sum 
% 
% Author: G. Richard, Janv. 2005 - MAJ:2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;


%Lecture du signal
[x,Fs]=audioread('sons_multipitch/A4_piano.wav');
%[x,Fs]=audioread('sons_multipitch/E4_oboe.wav'); 
%[x,Fs]=audioread('sons_multipitch/A3A4_piano.wav'); 
%[x,Fs]=audioread('sons_multipitch/A3C4E4G4_piano.wav'); 
soundsc(x,Fs);

spect_smooth=2; % options: 0: no spectral smoothness
               %          1: spectral smoothness principle with triangle 
               %           2 spectral smoothness principle with mean of three harmonics 

N=floor(0.7*Fs);      % Window size of analysed signal (only one window of signal is analysed)
Fmin=100;             % Minimal F0 frequency that can be detected
Fmax=900;             % Maximal F0 frequency that can be detected
H=4;                  % H = nombre de versions compress�es
prod_spec = 0;        %  m�thode pour la d�tection de pitch est produit spectral sinon par autocorrelation
freq_fond=[];         % tableau contenant la valeur des fr�quences fondamentales


%Minimal frequency resolution
dF_min=Fs/N;             

%beta coef d'inharmonicit�
beta=0.005; % OK A4 PIANO 
%alpha: coefficient for harmonic location around the true theoric harmony
alpha= 0.1;

%Window
w=hamming(N);

%Beginning of signal (e.g. attack) is discarded
offset=floor(0.1*Fs);
xw=x(offset+1:offset+N).*w;    %xw est la fenetre de signal analys�

%Minimal number of data points to satisfy the minimal frequency resolution
Nfft_min=Fs/dF_min;

%compute the smallest power of two that satisfies the minimal frequency resolution for FFT
p = nextpow2(Nfft_min);
Nfft=2^p;

%calcul FFT
Xk=fft(xw,Nfft);

%frequency resolution of FFT
df=Fs/Nfft;

%normalisation
Xk=Xk/max(abs(Xk)+eps);

%"Reduced" frequency
f=[0:Nfft-1]/Nfft+eps;

% fr�quence maximale
Rmax = floor((Nfft-1)/(2*H));



        figure(1);
        plot([1:Nfft/2]*(Fs/Nfft), 20*log10(abs(Xk(1:Nfft/2))),'b');
        axis([0 5000 -80 0]);
        xlabel('Fr�quences (Hz)','fontsize',14);
        ylabel('Spectre d''Amplitude (dB)','fontsize',14);      
        title('Spectre d''amplitude du signal original','fontsize',14);
        %pause

%loop on the number of pitches
%example criterion could use an energy ratio
i = 1;
seuil_F0 = 0.4;
criterion = 100;
while criterion > seuil_F0 
    
    %detection of main F0
        %Compute spectral sum
        % locate maximum
        %store value of estimated F0
    if(prod_spec == 1)
        S = ones(1,Rmax);
        for i=1:H
            S = S.*abs(Xk(1:i:i*Rmax));
        end
    else       
        S = zeros(1,Rmax);
        for i=1:H
            S = S + abs(Xk(1:i:i*Rmax));
        end
    end
    
    
    
    max_value = -1;
    F0 = -1;
    
    Nmin = floor(Fmin/df);
    Nmax = floor(Fmax/df);
    
    for i = Nmin:Nmax
        if S(i)>max_value
            max_value = S(i);
            F0 = df*(i-1);
        end
    end
    
    disp(F0);
    %Subtraction of main note (Main F0 with its harmonics)
        
        %localisation of harmonics around theoretical values (with or without inharmonicy coeeficient) 
        % beta: harmonicity coefficient ;  alpha: coefficient of tolerance
                
        % Harmonic suppression (wideness of an harmonic to be suppressed
        % depends on the main lob of the TF of the analysis window)
            % suppression of harmonics is done on abs(Xk) on forcing all values
            % of a harmonic peak to the minimum value of the peak (e.g. the
            % level of noise).
    k=1;
    while(k*F0/df < Rmax)
        max_value = -1;
        max_freq_harm = -1;
        f_inharm = k*F0*sqrt(1+k*k*beta);
        Nmin = floor((1-alpha)*f_inharm/df);
        Nmax = min(Nfft/2,floor((1+alpha)*f_inharm/df));
        
        for i = Nmin:Nmax
            if abs(Xk(i))> max_value
                max_value = abs(Xk(i));
                max_freq_harm = df*(i-1);
                max_index = i;
            end
        end

        k1 = max_index - floor(10/df) ;
        k2 = min(Nfft/2,max_index + floor(10/df)); 

        ang = angle(Xk(k1:k2));
        X = min(abs(Xk(k1:k2)));
        
        Xk(k1:k2) =  X*exp(1j*ang);
        k = k + 1;
    end
    
    criterion = max(abs(Xk(1:(Nfft/2))));
    
    figure(i);
    plot([1:Nfft/2]*(Fs/Nfft), 20*log10(abs(Xk(1:Nfft/2))),'b');
    axis([0 5000 -80 0]);
    xlabel('Frequences (Hz)','fontsize',14);
    ylabel('Spectre d''Amplitude (dB)','fontsize',14);      
    title('Spectre d''amplitude du signal original after harmonic suppression','fontsize',14);
    i = i+1; 
         
end

%end of loop