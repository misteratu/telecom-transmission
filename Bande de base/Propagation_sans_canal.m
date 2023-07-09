close all;
clear all;

Fe = 24000;         % Fréquence d'échantillonnage
Rb = 3000;          % Débit binaire
N_bits = 16;        % Nombre de bits
Tb = 1/Rb;          % Période de transmission d'un bit
Te = 1/Fe;          % Période d'échantillonnage
bits = randi([0, 1], 1 , N_bits);
Ns = fix(Tb/Te);    % Facteur de suréchantillonnage
NRZ = repelem(bits, 1, Ns);
Temps = linspace(0, N_bits/Rb, N_bits*Ns);
n0 = 8; % Taux d'erreur binaire à 0 lorsque n0 = 8 
        % Taux d'erreur binaire à 51 lorsque n0 = 3

%% Modulateur bande de bas avec mapping binaire à moyenne nulle

M = 2;
Rs1 = Rb / log2(M);
Ns = fix(Fe/Rs1);
h = ones(1, Ns);
I = find(bits == 0);
donnee = bits;

% Mapping :
% 1 -> 1
% 0 -> -1

donnee(I) = -1;
donnee = kron(donnee, [1 zeros(1, Ns - 1)]);

%% Filtre de mise en forme

x = filter(h, 1, donnee);

%% Passage par le canal 

bruit = 0; % On choisit de ne pas mettre de bruit
x_bruit = x + bruit;
hc = zeros(1, Ns);
hc(1) = 1;

x_canal = filter(hc, 1, x_bruit);

%% Filtre de réception
hr = h;
z = filter(hr, 1, x_canal);

figure;
plot(Temps, z);
xlabel("Temps en secondes (s)");
ylabel("signal en sortie du filtre de réception");

%% Echantilloneur : Démoduleur bande de base 
Mat = reshape(z, Ns, length(z)/Ns);
reception =  Mat(n0,:)/8; 
reception(reception <= 0) = 0;
reception(reception > 0) = 1;

% Calcul du taux d'erreur binaire
erreur = (reception == bits);
Taux_erreur = 100-100*mean(erreur)

% Calcul de g
g = conv(h, h);
g = conv(g, hc);

% Tracé de la réponse impulsionnelle globale de la chaine de transmission g
figure;
echelleg = (0:length(g)-1)*Tb/(length(g)-1);
plot(echelleg, g);
xlabel("Temps en secondes (s)");
ylabel("g(t)");

% tracé du diagramme de l'oeil en sortie du filtre de réception
figure;
plot(Mat);






