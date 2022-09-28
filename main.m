clear;
close all;
clc;

%% Tache 1 ----------------------------------------------------------------
Nb = 100000;
Ts = 1e-6; %temps symbole
fe = 20e6; %fréquence d'échantillonage
fse = 20;  %fréquence de suréchantillonage
Te = 1/fe; %période d'échantillonage
%%
b = randi([0 1], 1, Nb);
b_bis = b;
b = modulatePPM(b, fse);
b_bis = modulatePPM(b_bis, fse);
b = awgn(b, 10000);
%% Périodogramme de Welch
figure;
Nfft = 512;

%zéro padding
if Nfft>length(b)
    b = [b zeros(1, Nfft-length(b))];
end
    
fmin = -fe/2;
fmax = fe/2;
[f_axis, DSP] = Mon_Welch(b, Nfft, fe, fmin, fmax);
P1 = trapz(DSP);
% semilogy(f_axis, DSP);
%xlim([fmin, fmax]);

%% Spectre de puissance
%hold on;
TF = fftshift(fft(b, Nfft));
SP = abs(TF).^2/Nfft;
P2 = trapz(SP);
% semilogy(f_axis, SP);


%% Periodogramme de Daniel
%hold on;
window = 10;
DSP_daniel = Daniel(SP, window);
P3 = trapz(DSP_daniel);
% semilogy(f_axis, DSP_daniel);


% 

% figure;
% hold on;
% plot(f, P1*ones(1, Nfft));
% plot(f, P2*ones(1, Nfft));
% plot(f, P3*ones(1, Nfft));
% xlabel("Fréquence (Hz)");
% ylabel("Puissance");
% title("Puissance estimée par la méthode des trapèzes");

%% Méthode de Capon

DSP_th = @(f)  (((Ts^3)*((pi*f).^2))/16).*((sinc((Ts/2)*f)).^4);
semilogy(f_axis, 5*fe*DSP_th(f_axis));
ylim([10e-10, 10e5]);

Nfft = 1024;
[f_axis, DSP] = Mon_Welch(b, Nfft, fe, fmin, fmax);
TF = fftshift(fft(b, Nfft));
Rx = autocorr(b, Nfft-1); %matrice d'autocorrélation de b
Rx = toeplitz(Rx);
invRx = inv(Rx);

a = zeros(1, Nfft);
mat_a = zeros(Nfft);
index = 1;
for i = f_axis %on parcourt toutes les fréquences de f
    for n=1:Nfft
        a(n) = exp(-j*2*pi*(n-1)*i/fe);
    end
        mat_a(index, :) = a; %on stocke tous les vecteurs a dans une matrice
        index = index+1;
end

h_opti = zeros(1, Nfft);
mat_h_opti = zeros(Nfft);
index = 1;
for f=1:Nfft
    f
    h_opti = ((invRx*mat_a(f, :)')/(mat_a(f, :)*invRx*(mat_a(f, :)')))'; %on calcul le filtre optimal pour chaque fréquence
    mat_h_opti(f, :) = h_opti; %on stocke les filtres optimaux dans une matrice
end

SP1 = zeros(1, Nfft);
for f=1:Nfft
    SP1(f) = Nfft*abs(mat_h_opti(f, :)*Rx*mat_h_opti(f, :)'); %on calcule la puissance instantanée pour chaque fréquence
end

P4 = trapz(SP1);

hold on;
semilogy(f_axis, SP1);
xlabel("f (Hz)");
ylabel("Puissance");
title("DSP");
legend("Courbe théorique", "Méthode de Capon (lag=512)");
