close all;
clear all;

Fe = 24000;
Rb = 3000;
N_bits = 10000;
Tb = 1/Rb;
Te = 1/Fe;
bits = randi([0, 1], 1 , N_bits);
Ns = fix(Tb/Te);    % Facteur de suréchantillonnage
NRZ = repelem(bits, 1, Ns);
Temps = linspace(0, N_bits/Rb, N_bits*Ns);
n0 = 8; 
BW = 8000;

% Choix du rapport signal à bruit par bit souhaité à l entrée du récepteur
ensemble_R = linspace(0.001, 8, 10);
ensemble_TEB_exp1 = zeros(1, length(ensemble_R));


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

N = 99; 
i = 1; % Indice de parcours
for R = ensemble_R 
% Mise en place d'un bruit
    Px = mean(abs(x).^2);
    sigma = sqrt(Px*Ns/(2*log2(M)*10^(R/10)));
    bruit = sigma * randn(1, length(x));
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
    reception =  Mat(n0,:); 
    reception(reception <= 0) = 0;
    reception(reception > 0) = 1;
    
    erreur = (reception == bits);
    Taux_erreur = 1-mean(erreur);
    ensemble_TEB_exp1(i) = Taux_erreur;
    i = i + 1;
end 

n0 = 4; 
ensemble_TEB_exp2 = zeros(1, length(ensemble_R));

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

N = 99; 
i = 1; % Indice de parcours
for R = ensemble_R 
% Mise en place d'un bruit
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
    hr = h(1:Ns/2);
    z = filter(hr, 1, x_canal);

%% Echantilloneur : Démoduleur bande de base 
    Mat = reshape(z, Ns, length(z)/Ns);
    reception =  Mat(n0,:); 
    reception(reception <= 0) = 0;
    reception(reception > 0) = 1;
    
    erreur = (reception == bits);
    Taux_erreur = 1-mean(erreur);
    ensemble_TEB_exp2(i) = Taux_erreur;
    i = i + 1;
end 


n0 = 16; 
ensemble_TEB_exp3 = zeros(1, length(ensemble_R));


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
    reception =  Mat(n0,:)/n0; 
    reception_demapee = zeros(1, 2*length(Mat(n0,:)));

    % Démodulation mapping classique 

%     I00 = find(reception <= -2);
%     I10 = find(reception <= 0 & reception > -2);
%     I11 = find(reception > 2);
%     I01 = find(reception <= 2 & reception > 0);
%     
%     reception_demapee(2*I01) = 1;
%     reception_demapee(2*I11-1) = 1;
%     reception_demapee(2*I11) = 1;
%     reception_demapee(2*I10-1) = 1;

    % Démodulation mapping de gray

    I00 = find(reception <= -2);
    I01 = find(reception <= 0 & reception > -2);
    I10 = find(reception > 2);
    I11 = find(reception <= 2 & reception > 0);
    
    reception_demapee(2*I01) = 1;
    reception_demapee(2*I11-1) = 1;
    reception_demapee(2*I11) = 1;
    reception_demapee(2*I10-1) = 1;
    
    erreur = (reception_demapee == bits);
    Taux_erreur = 1 - mean(erreur);
    ensemble_TEB_exp3(i) = Taux_erreur;
    i = i + 1;
end 

n0 = 16; 
ensemble_TEB_exp4 = zeros(1, length(ensemble_R));


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

%% Filtre de mise en forme

h = ones(1, Ns);
x = filter(h, 1, kron_pair);


%% Passage par le canal 

N = 99; 
i = 1; % Indice de parcours
for R = ensemble_R 
    
    % Mise en place d'un bruit
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
    reception =  Mat(n0,:)/n0; 
    reception_demapee = zeros(1, 2*length(Mat(n0,:)));

    % Démodulation mapping classique 

     I00 = find(reception <= -2);
     I10 = find(reception <= 0 & reception > -2);
     I11 = find(reception > 2);
     I01 = find(reception <= 2 & reception > 0);
     
     reception_demapee(2*I01) = 1;
     reception_demapee(2*I11-1) = 1;
     reception_demapee(2*I11) = 1;
     reception_demapee(2*I10-1) = 1;

    
    
    erreur = (reception_demapee == bits);
    Taux_erreur = 1 - mean(erreur);
    ensemble_TEB_exp4(i) = Taux_erreur;
    i = i + 1;
end 


%% Affichage des différentes comparaisons des TEB des 3 chaines.
figure;
semilogy(ensemble_R, ensemble_TEB_exp1,'LineWidth',2);
hold on
semilogy(ensemble_R, ensemble_TEB_exp2,'LineWidth',2);
hold off
xlabel("Rapport signal à bruit par bit souhaité : R");
ylabel("Taux d'erreur binaire");
legend("Taux d'erreur experimental chaine 1", "Taux d'erreur experimental chaine 2");

figure;
semilogy(ensemble_R, ensemble_TEB_exp1,'LineWidth',2);
hold on
semilogy(ensemble_R, ensemble_TEB_exp3,'LineWidth',2);
hold off
xlabel("Rapport signal à bruit par bit souhaité : R");
ylabel("Taux d'erreur binaire");
legend("Taux d'erreur experimental chaine 1", "Taux d'erreur experimental chaine 3");

ensemble_TEB_th = 2*(M-1)/M*qfunc(sqrt((6*log2(M))/(M*M-1)*10.^(ensemble_R/10)))/log2(M);

figure;
semilogy(ensemble_R, ensemble_TEB_exp3,'LineWidth',2);
hold on
semilogy(ensemble_R, ensemble_TEB_exp4,'LineWidth',2);
semilogy(ensemble_R, ensemble_TEB_th,'LineWidth',2);
hold off
xlabel("Rapport signal à bruit par bit souhaité : R");
ylabel("Taux d'erreur binaire");
legend("Taux d'erreur experimental chaine 3 Gray", "Taux d'erreur experimental chaine 3 Non Gray", "Taux d'erreur experimental chaine 3 théorique");

