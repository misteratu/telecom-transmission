close all;
clear all;

Fe = 24000;     % Fréquence d'échantillonnage
fp = 2000;      % Fréquence de la porteuse
Rb = 3000;      % Débit binaire
N_bits = 3000;   % Nombre de bits transmis
Tb = 1/Rb;      % Période de transmission d'un bit
Te = 1/Fe;      % Période d'échantillonnage 
Ns = fix(Tb/Te);    % Facteur de suréchantillonnage
n0 = 1;
M = 8;

% Choix du rapport signal à bruit par bit souhaité à l entrée du récepteur
ensemble_R = linspace(0.001, 6, 10);
ensemble_TEB_exp = zeros(1, length(ensemble_R));

bits = randi([0, 1], 1 , N_bits);   % bits d'information à transmettre 
L = 10;
bits_col = reshape(bits, 3, length(bits)/3);
var_map = bi2de(bits_col.');
m = pskmod(var_map, 8, 0, 'gray');
m = m.';


%% Créer la réponse impulsionnelle avec des symboles binaires à moyenne nulle.


Rs3 = Rb / log2(M);
Ts3 = 1/Rs3;
Ns3 = fix(Fe/Rs3);
alpha = 0.20;
h = rcosdesign(alpha, L, M, 'sqrt');
kron_cos = kron(m, [1 zeros(1,Ns3 - 1)]);

retard = M*L/2;
kron_cos = [kron_cos, kron_cos(1:retard)]; % Gestion du retard

signal_sortie = filter(h, 1, kron_cos);

x = signal_sortie(retard+1:length(signal_sortie));


%% Passage par le canal 

i = 1; % Indice de parcours
for R = ensemble_R 
    Px = mean(abs(x).^2);
    sigma = sqrt(Px*Ns3/(2*log2(M)*10^(R/10)));
    bruit_i = sigma * randn(1, length(x));
    bruit_q = sigma * randn(1, length(x));
    bruit = bruit_i + 1j*bruit_q;
    %bruit = 0;  % A modifier si on veut du bruit
    x_bruit = x + bruit;
    
    %% Démodulateur

    x_dem = [x_bruit, x_bruit(1:retard)]; % Gestion du retard
    signal_sortie = filter(h, 1, x_dem);
    signal_sortie = signal_sortie(retard+1:length(signal_sortie));
    
    Mat = reshape(signal_sortie, Ns3, length(signal_sortie)/Ns3);
 
    reception =  Mat(n0,:);
    
    reception_demod = pskdemod(reception, M, 0, 'gray');


    I0 = find(reception_demod == 0);
    I1 = find(reception_demod == 1);
    I2 = find(reception_demod == 2);
    I3 = find(reception_demod == 3);
    I4 = find(reception_demod == 4);
    I5 = find(reception_demod == 5);
    I6 = find(reception_demod == 6);
    I7 = find(reception_demod == 7);

    bits_reception = zeros(1, N_bits);
    bits_reception(3*I1-2) = 1;

    bits_reception(3*I2-1) = 1;

    bits_reception(3*I3-2) = 1;
    bits_reception(3*I3-1) = 1;

    bits_reception(3*I4) = 1;

    bits_reception(3*I5) = 1;
    bits_reception(3*I5-2) = 1;

    bits_reception(3*I6) = 1;
    bits_reception(3*I6-1) = 1;

    bits_reception(3*I7-2) = 1;
    bits_reception(3*I7-1) = 1;
    bits_reception(3*I7) = 1;

    erreur = (bits_reception == bits);
    % Calcul du taux d'erreur binaire pour une valeur de Eb/n0
    Taux_erreur = 1-mean(erreur);
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
