close all;
clear all;

Fe = 24000;     % Fréquence d'échantillonnage
Rb = 3000;      % Débit binaire
N_bits = 300;   % Nombre de bits transmis
Tb = 1/Rb;      % Période de transmission d'un bit
Te = 1/Fe;      % Période d'échantillonnage 
Ns = fix(Tb/Te);    % Facteur de suréchantillonnage

bits = randi([0, 1], 1 , N_bits);   % bits d'information à transmettre 
NRZ = repelem(bits, 1, Ns);     % Singal NRZ des bits

Temps = linspace(0, N_bits/Rb, N_bits*Ns);
L = 10;

%% Créer la réponse impulsionnelle avec des symboles binaires à moyenne nulle.

M = 2;
Rs1 = Rb / log2(M);
Ts1 = 1/Rs1;
Ns1 = fix(Fe/Rs1);
h = ones(1, Ns1);
I = find(bits == 0);
donnee1 = bits;

% Mapping :
% 1 -> 1
% 0 -> -1

donnee1(I) = -1;
donnee1 = kron(donnee1, [1 zeros(1, Ns1 - 1)]);
signal_sortie1 = filter(h, 1, donnee1);

% Calcul de la densité spectrale du signal de sortie
[DSP1, F1] = pwelch(signal_sortie1, [], [], [], Fe, 'twosided');
Freq1 = linspace(-Fe/2, Fe/2, length(F1));
Temps1 = linspace(0, N_bits/Rb, N_bits*Ns1);

%% Créer la réponse impulsionnelle avec des symboles 4-aires à moyenne nulle.

M = 4;
Rs2 = Rb / log2(M);
Ts2 = 1/Rs2;
Ns2 = fix(Fe/Rs2);

% Mapping :
% 11 -> 3
% 01 -> 1
% 10 -> -1
% 00 -> -3

donnee_pair  = (2*bi2de(reshape(bits, 2, length(bits)/2).')-3).';
kron_pair = kron(donnee_pair, [1 zeros(1, Ns2 - 1)]);
h = ones(1, Ns2);
signal_sortie2 = filter(h, 1, kron_pair);

% Calcul de la densité spectrale du signal de sortie
[DSP2, F2] = pwelch(signal_sortie2, [], [], [], Fe, 'twosided');
Freq2 = linspace(-Fe/2, Fe/2, length(F2));
Temps2 = linspace(0, N_bits/Rb, N_bits*Ns2/2);

%% Créer la réponse impulsionnelle avec des symboles binaires à moyenne nulle.
M = 2;
Rs3 = Rb / log2(M);
Ts3 = 1/Rs3;
Ns3 = fix(Fe/Rs3);
alpha = 0.5;
h = rcosdesign(alpha, L, Ns);
symboles=-2*bits-1;
kron_cos = kron(symboles, [1 zeros(1,Ns3 - 1)]);
signal_sortie3 = filter(h, 1, kron_cos);

% Calcul de la densité spectrale du signal de sortie
[DSP3, F3] = pwelch(signal_sortie3, [], [], [], Fe, 'twosided');
Freq3 = linspace(-Fe/2, Fe/2, length(F3));
Temps3 = linspace(0, N_bits/Rb, N_bits*Ns3);

%% Tracé des signaux générés en sortie des filtres

figure;
subplot(3,1,1);
plot(Temps1, signal_sortie1)
xlabel("Temps en secondes (s)");
ylabel("Signal modulé");
title("Tracé du signal modulé du modulateur 1");

subplot(3,1,2);
plot(Temps2, signal_sortie2)
xlabel("Temps en secondes (s)");
ylabel("Signal modulé");
title("Tracé du signal modulé du modulateur 2");

subplot(3,1,3);
plot(Temps3, signal_sortie3)
xlabel("Temps en secondes (s)");
ylabel("Signal modulé");
title("Tracé du signal modulé du modulateur 3");

%% Densité spectrale théorique du modulateur 1
Sx1 = Ts1*(sinc(Freq1*Ts1)).^2;

%% Densité spectrale théorique du modulateur 2
Sx2 = Ts2*(sinc(Freq2*Ts2)).^2;

%% Densité spectrale théorique du modulateur 3
I1 = find(abs(Freq3) <= (1-alpha)/(2*Ts3));
I2 = find(abs(Freq3) >= (1-alpha)/(2*Ts3) & abs(Freq3) <= (1+alpha)/(2*Ts3));
Sx3 = zeros(1, length(Freq3));
Sx3(I1) = Ts3;
theta = (pi*Ts3)/alpha * (abs(Freq3(I2)) - (1-alpha)/(2*Ts3));
Sx3(I2) = (Ts3/2)*(1+cos(theta));

%% Tracés des densités spectrales numériques 
figure;
subplot(3,1,1);
semilogy(Freq1,fftshift(DSP1));
xlabel("Fréquence en hertz (Hz)");
ylabel("DSP du signal");
title("Tracé de la densité spectrale du signal modulé du modulateur 1");

subplot(3,1,2);
semilogy(Freq2,fftshift(DSP2));
xlabel("Fréquence en hertz (Hz)");
ylabel("DSP du signal");
title("Tracé de la densité spectrale du signal modulé du modulateur 2");

subplot(3,1,3);
semilogy(Freq3,fftshift(DSP3));
xlabel("Fréquence en hertz (Hz)");
ylabel("DSP du signal");
title("Tracé de la densité spectrale du signal modulé du modulateur 3");


%% Tracés des densités spectrales théoriques et numériques 
figure;
subplot(3,1,1);
semilogy(Freq1,fftshift(DSP1));
hold on 
semilogy(Freq1, Sx1);
hold off
xlabel("Fréquence en hertz (Hz)");
ylabel("DSP du signal");
title("Tracé de la densité spectrale théorique et numérique du signal modulé du modulateur 1");
legend("DSP numérique", "DSP théorique");

subplot(3,1,2);
semilogy(Freq2,fftshift(DSP2));
hold on 
semilogy(Freq2, Sx2);
hold off
xlabel("Fréquence en hertz (Hz)");
ylabel("DSP du signal");
title("Tracé de la densité spectrale théorique et numérique du signal modulé du modulateur 2");
legend("DSP numérique", "DSP théorique");

subplot(3,1,3);
semilogy(Freq3,fftshift(DSP3));
hold on
semilogy(Freq3, Sx3);
hold off
xlabel("Fréquence en hertz (Hz)");
ylabel("DSP du signal");
title("Tracé de la densité spectrale théorique et numérique du signal modulé du modulateur 3");
legend("DSP numérique", "DSP théorique");


%% Tracés de la superposition de la densité spectrale du signal des 3 modulateurs

figure;
semilogy(Freq1,fftshift(DSP1));
hold on
semilogy(Freq2,fftshift(DSP2));
semilogy(Freq3,fftshift(DSP3));
hold off
xlabel("Fréquence en hertz (Hz)");
ylabel("DSP des signaux en sortie");
legend("DSP modulateur 1", "DSP modulateur 2", "DSP modulateur 3");