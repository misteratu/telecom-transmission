close all;
clear all;

Fe = 24000;         % Fréquence d'échantillonnage
Rb = 3000;          % Débit binaire
N_bits = 10000;     % Nombre de bits transmis
Tb = 1/Rb;          % Période de transmission d'un bit
Te = 1/Fe;          % Période d'échantillonnage 
bits = randi([0, 1], 1 , N_bits);   % bits d'information à transmettre
Ns = fix(Tb/Te);    % Facteur de suréchantillonnage
NRZ = repelem(bits, 1, Ns);
Temps = linspace(0, N_bits/Rb, N_bits*Ns);
n0 = 16; % Taux d'erreur binaire à 0 lorsque n0 = 8 
        % Taux d'erreur binaire à 51 lorsque n0 = 3
BW = 8000;

% Choix du rapport signal à bruit par bit souhaité à l entrée du récepteur
ensemble_R = linspace(0.001, 8, 10);
ensemble_TEB_exp = zeros(1, length(ensemble_R));
% R = Eb/n0 


%% Modulateur bande de bas avec mapping binaire à moyenne nulle

M = 4;
Rs2 = Rb / log2(M);
Ts2 = 1/Rs2;
Ns = fix(Fe/Rs2);

% Mapping :
% 11 -> 3
% 01 -> 1
% 10 -> -1
% 00 -> -3

donnee_pair  = (2*bi2de(reshape(bits, 2, length(bits)/2).')-3).';
kron_pair = kron(donnee_pair, [1 zeros(1, Ns - 1)]);

% Modification du mapping pour qu'il soit de gray

Igray10 = find(kron_pair == -1);
Igray11 = find(kron_pair == 3);
Igray01 = find(kron_pair == 1);

kron_pair(Igray01) = -1;
kron_pair(Igray11) = 1;
kron_pair(Igray10) = 3;

%% Filtre de mise en forme

h = ones(1, Ns);
x = filter(h, 1, kron_pair);


%% Passage par le canal 

N = 99; 
i = 1; % Indice de parcours
for R = ensemble_R 
    
    % Mise en place d'un bruit
    %R = 0;
    Px = mean(abs(x).^2);
    sigma = sqrt(Px*Ns/(2*log2(M)*10^(R/10)));
    bruit = sigma * randn(1, length(x));
    %bruit = 0;
    x_bruit = x + bruit;
    
    % Filtre du filtre  
    x_bruit = [x_bruit, x_bruit(1:N)]; % Gestion du retard en entrée du canal
    Taille_filtre = -N:1:N;
    Passe_bas_i = 2*BW/Fe*sinc(2*BW/Fe*Taille_filtre);
    
    % Gestion du retard en sortie du canal
    x_canal = filter(Passe_bas_i, 1, x_bruit);
    x_canal = x_canal(N+1:length(x_canal));
    
    %% Filtre de réception
    hr = h;
    z = filter(hr, 1, x_canal);
    
    %% Echantilloneur : Démoduleur bande de base 

    Mat = reshape(z, Ns, length(z)/Ns);
    
    figure;
    plot(Mat);
    title("Diagramme de l'oeil pour R = ", num2str(R));

    reception =  Mat(n0,:)/n0; 
    reception_demapee = zeros(1, 2*length(Mat(n0,:)));

    %% Démodulation mapping classique 

%     I00 = find(reception <= -2);
%     I10 = find(reception <= 0 & reception > -2);
%     I11 = find(reception > 2);
%     I01 = find(reception <= 2 & reception > 0);
%     
%     reception_demapee(2*I01) = 1;
%     reception_demapee(2*I11-1) = 1;
%     reception_demapee(2*I11) = 1;
%     reception_demapee(2*I10-1) = 1;

    %% Démodulation mapping de gray

    I00 = find(reception <= -2);
    I01 = find(reception <= 0 & reception > -2);
    I10 = find(reception > 2);
    I11 = find(reception <= 2 & reception > 0);
    
    reception_demapee(2*I01) = 1;
    reception_demapee(2*I11-1) = 1;
    reception_demapee(2*I11) = 1;
    reception_demapee(2*I10-1) = 1;
    
    % Calcul du taux d'erreur binaire pour une valeur de Eb/n0
    erreur = (reception_demapee == bits);
    Taux_erreur = 1 - mean(erreur);
    ensemble_TEB_exp(i) = Taux_erreur;
    i = i + 1;
end 

% Calcul du taux d'erreur binaire théorique
ensemble_TEB_th = 2*(M-1)/M*qfunc(sqrt((6*log2(M))/(M*M-1)*10.^(ensemble_R/10)))/log2(M);

%% Affichage des courbes de taux d'erreurs binaire
figure;
semilogy(ensemble_R, ensemble_TEB_exp,'LineWidth',2);
hold on
semilogy(ensemble_R, ensemble_TEB_th,'LineWidth',2);
hold off
xlabel("Rapport signal à bruit par bit souhaité : R");
ylabel("Taux d'erreur binaire");
legend("Taux d'erreur experimental", "Taux d'erreur théorique");
