close all;
clear all;

Fe = 24000;         % Fréquence d'échantillonnage
Rb = 6000;          % Débit binaire
N_bits = 30000;     % Nombre de bits transmis
Tb = 1/Rb;          % Période de transmission d'un bit
Te = 1/Fe;          % Période d'échantillonnage 
bits = randi([0, 1], 1 , N_bits);   % bits d'information à transmettre 
Ns = fix(Tb/Te);    % Facteur de suréchantillonnage
NRZ = repelem(bits, 1, Ns);
Temps = linspace(0, N_bits/Rb, N_bits*Ns);
n0 = 4; % Taux d'erreur binaire à 0 lorsque n0 = 8 
        % Taux d'erreur binaire à 51 lorsque n0 = 3

phi = deg2rad(100);



%% Sans correction Erreur phase

% Choix du rapport signal à bruit par bit souhaité à l entrée du récepteur
ensemble_R = linspace(0.001, 6, 10);
ensemble_TEB_expSEP = zeros(1, length(ensemble_R));


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

donneeCode = donnee;
for r = 2:N_bits
    donneeCode(r) = donnee(r)*donneeCode(r-1);
end

donnee = kron(donneeCode, [1 zeros(1, Ns - 1)]);

%% Filtre de mise en forme

x = filter(h, 1, donnee);


%% Passage par le canal 

N = 99; 
i = 1; % Indice de parcours
for R = ensemble_R 
% Mise en place d'un bruit
    Px = mean(abs(x).^2);
    sigma = sqrt(Px*Ns/(2*log2(M)*10^(R/10)));
    
    bruit_i = sigma * randn(1, length(x));
    bruit_q = sigma * randn(1, length(x));
    bruit = bruit_i + 1j*bruit_q;

    %bruit = 0;  % A modifier si on veut du bruit
    x_bruit = x + bruit;

    x_bruit = x_bruit.*exp(1j*phi);

%% Filtre de réception
    hr = h;
    z = filter(hr, 1, x_bruit);

%% Echantilloneur : Démoduleur bande de base 
    Mat = reshape(z, Ns, length(z)/Ns);
    reception =  Mat(n0,:);

    reception = sign(real(reception));

    receptionCode = reception;
    for r = 2:N_bits
        receptionCode(r) = reception(r)*reception(r-1);
    end

    receptionCode(receptionCode <= 0) = 0;
    receptionCode(receptionCode > 0) = 1;
    
    % Calcul du taux d'erreur binaire pour une valeur de Eb/n0
    erreur = (receptionCode == bits);
    Taux_erreur = 1-mean(erreur);
    ensemble_TEB_expSEP(i) = Taux_erreur;
    i = i + 1;
end 

%% Sans codage
% Choix du rapport signal à bruit par bit souhaité à l entrée du récepteur
ensemble_R = linspace(0.001, 6, 10);
ensemble_TEB_expSC = zeros(1, length(ensemble_R));


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

donneeCode = donnee;

donnee = kron(donneeCode, [1 zeros(1, Ns - 1)]);

%% Filtre de mise en forme

x = filter(h, 1, donnee);


%% Passage par le canal 

N = 99; 
i = 1; % Indice de parcours
for R = ensemble_R 
% Mise en place d'un bruit
    Px = mean(abs(x).^2);
    sigma = sqrt(Px*Ns/(2*log2(M)*10^(R/10)));
    
    bruit_i = sigma * randn(1, length(x));
    bruit_q = sigma * randn(1, length(x));
    bruit = bruit_i + 1j*bruit_q;

    %bruit = 0;  % A modifier si on veut du bruit
    x_bruit = x + bruit;

    x_bruit = x_bruit.*exp(1j*phi);

%% Filtre de réception
    hr = h;
    z = filter(hr, 1, x_bruit);

%% Echantilloneur : Démoduleur bande de base 
    Mat = reshape(z, Ns, length(z)/Ns);
    reception =  Mat(n0,:);

    phi_cor = 0.5*angle(sum(reception.^2));
    reception = reception*exp(-1j*phi_cor);

    reception = sign(real(reception));

    receptionCode = reception;

    receptionCode(receptionCode <= 0) = 0;
    receptionCode(receptionCode > 0) = 1;
    
    % Calcul du taux d'erreur binaire pour une valeur de Eb/n0
    erreur = (receptionCode == bits);
    Taux_erreur = 1-mean(erreur);
    ensemble_TEB_expSC(i) = Taux_erreur;
    i = i + 1;
end 


%% Avec les deux 

% Choix du rapport signal à bruit par bit souhaité à l entrée du récepteur
ensemble_R = linspace(0.001, 6, 10);
ensemble_TEB_exp = zeros(1, length(ensemble_R));


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

donneeCode = donnee;
for r = 2:N_bits
    donneeCode(r) = donnee(r)*donneeCode(r-1);
end

donnee = kron(donneeCode, [1 zeros(1, Ns - 1)]);

%% Filtre de mise en forme

x = filter(h, 1, donnee);


%% Passage par le canal 

N = 99; 
i = 1; % Indice de parcours
for R = ensemble_R 
% Mise en place d'un bruit
    Px = mean(abs(x).^2);
    sigma = sqrt(Px*Ns/(2*log2(M)*10^(R/10)));
    
    bruit_i = sigma * randn(1, length(x));
    bruit_q = sigma * randn(1, length(x));
    bruit = bruit_i + 1j*bruit_q;

    %bruit = 0;  % A modifier si on veut du bruit
    x_bruit = x + bruit;

    x_bruit = x_bruit.*exp(1j*phi);

%% Filtre de réception
    hr = h;
    z = filter(hr, 1, x_bruit);

%% Echantilloneur : Démoduleur bande de base 
    Mat = reshape(z, Ns, length(z)/Ns);
    reception =  Mat(n0,:);

    phi_cor = 0.5*angle(sum(reception.^2));
    reception = reception*exp(-1j*phi_cor);

    reception = sign(real(reception));

    receptionCode = reception;
    for r = 2:N_bits
        receptionCode(r) = reception(r)*reception(r-1);
    end

    receptionCode(receptionCode <= 0) = 0;
    receptionCode(receptionCode > 0) = 1;
    
    % Calcul du taux d'erreur binaire pour une valeur de Eb/n0
    erreur = (receptionCode == bits);
    Taux_erreur = 1-mean(erreur);
    ensemble_TEB_exp(i) = Taux_erreur;
    i = i + 1;
end 


%% Affichage des courbes de taux d'erreurs binaire
figure;
semilogy(ensemble_R, ensemble_TEB_exp,'LineWidth',2);
hold on
semilogy(ensemble_R, ensemble_TEB_expSEP,'LineWidth',2);
semilogy(ensemble_R, ensemble_TEB_expSC,'LineWidth',2);
hold off
xlabel("Rapport signal à bruit par bit souhaité : R");
ylabel("Taux d'erreur binaire");
legend("TEB experimental avec erreur de phase de 100° corrigée et codage", "TEB experimental avec erreur de phase de 100° non corrigée et codage", "TEB experimental avec erreur de phase de 100° corrigée et sans codage");

