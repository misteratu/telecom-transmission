close all;
clear all;

Fe = 24000;     % Fréquence d'échantillonnage
fp = 2000;      % Fréquence de la porteuse
Rb = 3000;      % Débit binaire
N_bits = 10000;   % Nombre de bits transmis
Tb = 1/Rb;      % Période de transmission d'un bit
Te = 1/Fe;      % Période d'échantillonnage 
Ns = fix(Tb/Te);    % Facteur de suréchantillonnage
n0 = 1;
% Choix du rapport signal à bruit par bit souhaité à l entrée du récepteur
ensemble_R = linspace(0.001, 6, 10);
ensemble_TEB_exp = zeros(1, length(ensemble_R));

bits = randi([0, 1], 1 , N_bits);   % bits d'information à transmettre 
L = 10;
bits2 = bits;
bits2(bits == 0) = -1;
m1 = bits2(1:2:N_bits);
m2 = bits2(2:2:N_bits);

figure;
subplot(2,1,1);
plot(m1);
subplot(2,1,2);
plot(m2);

m = m1 + m2 * 1i;

%% Créer la réponse impulsionnelle avec des symboles binaires à moyenne nulle.

M = 4;
Rs3 = Rb / log2(M);
Ts3 = 1/Rs3;
Ns3 = fix(Fe/Rs3);
alpha = 0.35;
h = rcosdesign(alpha, L, Ns3, 'sqrt');
kron_cos = kron(m, [1 zeros(1,Ns3 - 1)]);

retard = Ns3*L/2;
kron_cos = [kron_cos, kron_cos(1:retard)]; % Gestion du retard

signal_sortie = filter(h, 1, kron_cos);

x = signal_sortie(retard+1:length(signal_sortie));
t = linspace(0, N_bits/Rb, N_bits*Ns3/2);

%% Calcul et affichage des DSP en quadrature et en phase

% Calcul
m1bis = filter(h, 1, m1);
m2bis = filter(h, 1, m2);

[DSP1, F1] = pwelch(m1bis, [], [], [], Fe, 'twosided');
[DSP2, F2] = pwelch(m2bis, [], [], [], Fe, 'twosided');
Freq1 = linspace(-Fe/2, Fe/2, length(F1));
Freq2 = linspace(-Fe/2, Fe/2, length(F2));

% Affichage
figure;
subplot(2,1,1);
semilogy(Freq1,fftshift(DSP1),'LineWidth',2);
xlabel("Fréquence en hertz (Hz)");
ylabel("DSP du signal");
title("Tracé de la densité spectrale du signal sur la voie en phase");

subplot(2,1,2);
semilogy(Freq2,fftshift(DSP2),'LineWidth',2);
xlabel("Fréquence en hertz (Hz)");
ylabel("DSP du signal");
title("Tracé de la densité spectrale du signal sur la voie en quadrature");

%% Affichage du signal en sortie x et de sa DSP

figure;
plot(t, x);

[DSPx, Fx] = pwelch(x, [], [], [], Fe, 'twosided');
Freqx = linspace(-Fe/2, Fe/2, length(Fx));

figure;
semilogy(Freqx,fftshift(DSPx),'LineWidth',2);
xlabel("Fréquence en hertz (Hz)");
ylabel("DSP du signal");
title("Tracé de la densité spectrale du signal transmis sur la fréquence porteuse");



%% Passage par le canal 

i = 1; % Indice de parcours
for R = ensemble_R 
    Px = mean(abs(x).^2);
    sigma = sqrt(Px*Ns3/(2*log2(M)*10^(R/10)));
    bruit_i = sigma * randn(1, length(x));
    bruit_q = sigma * randn(1, length(x));
    bruit = bruit_i + 1j*bruit_q;
    %bruit = 0;  % A modifier si on veut du bruit
    disp(sigma)
    x_bruit = x + bruit;
    
    %% Démodulateur

    x_dem = [x_bruit, x_bruit(1:retard)]; % Gestion du retard
    signal_sortie = filter(h, 1, x_dem);
    signal_sortie = signal_sortie(retard+1:length(signal_sortie));
    
    Mat = reshape(signal_sortie, Ns3, length(signal_sortie)/Ns3);
 
    reception =  Mat(n0,:); 
    
    reel = real(reception);
    im = imag(reception);
    
    %On vérifie le signe de reel ou de im pour savoir si c'est 0 ou 1
    
    im(im > 0) = 1; 
    im(im <= 0) = 0;
    
    reel(reel > 0) = 1;
    reel(reel <= 0) = 0;
    
    bits_reception = zeros(1, N_bits);
    bits_reception(1:2:N_bits) = reel;
    bits_reception(2:2:N_bits) = im;

    erreur = (bits_reception == bits);
    % Calcul du taux d'erreur binaire pour une valeur de Eb/n0
    Taux_erreur = 1-mean(erreur)
    ensemble_TEB_exp(i) = Taux_erreur;
    i = i + 1;
end 

figure;
plot(Mat);

scatterplot(reception);
TEB = 2*qfunc(sqrt(2*log2(M)*10.^(ensemble_R/10))*sin(pi/M))/log2(M);

figure;
semilogy(ensemble_R, ensemble_TEB_exp,'LineWidth',2);
hold on
semilogy(ensemble_R, TEB,'LineWidth',2);
hold off
