clear all ;
close all ;
clc;

%% Lecture

img = imread('code_barre2.png');
[h, width, c] = size(img);
fprintf("L'image est une matrice de %s entre %d et %d\n", class(img(10,10,2)), min(img(:)), max(img(:)));
fprintf("Elle est de dimension: %d*%d*%d\n", h, width, c);

figure, imshow(img);
[x, y] = ginput(2);
fprintf('les points choisies sont (%i, %i) et (%i, %i)\n', x(1), y(1), x(2), y(2));
% x = [443;69];
% y = [45;390];

%% Lancer Aléatoire d'un rayon
Ln = sqrt( (x(1)-x(2))^2 + (y(1)-y(2))^2 ); % longueur du segment

nbr_pt = 2*ceil(Ln); % Nombre de point a prendre en respectant le Crtr2Shannon

Mu = zeros(nbr_pt, 2);
Mu(1,:) = [x(1) y(1)];
Mu(nbr_pt,:) = [x(2) y(2)];

signature_img(1) = sum( img( ceil(Mu(1,2)), ceil(Mu(1,1)), :) )/3;

for u=2:nbr_pt-1
    Mu(u,:) = Mu(1,:) + (u/nbr_pt)*(Mu(nbr_pt,:) - Mu(1,:));
    signature_img(u) = sum( img( ceil(Mu(u,2)), ceil(Mu(u,1)), :) )/3;
end
signature_img(nbr_pt) = sum( img( ceil(Mu(nbr_pt,2)), ceil(Mu(nbr_pt,1)), :) )/3;
signature_img = signature_img(:);

% figure, plot(signature_img);
% title("Signature de l'image");

%% Extraction de la signature le long d’un rayon

max_signature = max(signature_img);
min_signature = min(signature_img);
interval_signature = min_signature:(max_signature - min_signature)/255:max_signature;


% figure, hist(signature_img, 256);
% title("Histogramme des valeurs de la signature");
hist_img = hist(signature_img, 256);

N = 256;
w = zeros(1,N);
mu = zeros(1,N);
denum = sum(hist_img);

for i=1:N
    w(i) = (sum(hist_img(1:i))) / denum;
    mu(i) = (sum([1:i] .* hist_img(1:i))) / denum;
end

crit = w.*((mu(1,N)-mu).^2) + (1 - w).*mu.^2;


[crit_max, k_max] = max(crit);


figure, plot(signature_img);
hold on, plot(interval_signature(k_max)*ones(nbr_pt));
title("Signature de l'image avec le seuil");

signature_binaire = signature_img > k_max;
figure, plot(signature_binaire);
ylim([-0.25 1.25]);
title("Signature Binaire de l'image");

% Estimation des points limites virtuels correspondants au premier 
% et au dernier échantillon de la signature utile

for i=1:nbr_pt
    if(signature_binaire(i) == 0)
        indice_premier = i;
        break;
    end
end
for i=0:nbr_pt-1
    if(signature_binaire(nbr_pt - i) == 0)
        indice_dernier = nbr_pt - i;
        break;
    end
end
fprintf("Indice du premier échantillon: %d, indice du dernier: %d\n", indice_premier, indice_dernier);

[x1, y1] = deal(Mu(indice_premier, 1), Mu(indice_premier, 2));
[x2, y2] = deal(Mu(indice_dernier, 1), Mu(indice_dernier, 1));

%% Nouvelle extraction de la signature binaire, cette fois entre le premier et dernier echantillons utiles

Ln2 = sqrt( (x1-x2)^2 + (y1-y2)^2 ); % longueur du segment

nbr_pt2 = round(2*ceil(Ln2)/95)*95;  % Nbr de pt en respectant le Crtr2Shannon et L=95u

Mu2 = zeros(nbr_pt2, 2);
Mu2(1,:) = [x1 y1];
Mu2(nbr_pt2,:) = [x2 y2];

signature_img2(1) = sum( img( ceil(Mu2(1,2)), ceil(Mu2(1,1)), :) )/3;

for u=2:nbr_pt2-1
    Mu2(u,:) = Mu2(1,:) + (u/nbr_pt2)*(Mu2(nbr_pt2,:) - Mu2(1,:));
    signature_img2(u) = sum( img( ceil(Mu2(u,2)), ceil(Mu2(u,1)), :) )/3;
end
signature_img2(nbr_pt2) = sum( img( ceil(Mu(nbr_pt2,2)), ceil(Mu(nbr_pt2,1)), :) )/3;
signature_img2 = signature_img2(:);

figure, plot(signature_img2);
hold on, plot(interval_signature(k_max)*ones(nbr_pt2));
title("Signature de l'image entre les points utiles avec le seuil");

signature_binaire2 = signature_img2 > k_max;
figure, plot(signature_binaire2);
ylim([-0.25 1.25]);
title("Signature Binaire de l'image entre les points utiles");

%% Identification des chiffres
unite_base_u = nbr_pt2/95; % l’unité de base u

% Décomposer la signature en ses différentes parties
garde_gauche   = signature_binaire2(1:3*unite_base_u);
region_1       = signature_binaire2(4*unite_base_u:45*unite_base_u); % représente les chiffres de rang 2 à 7
garde_centrale = signature_binaire2(46*unite_base_u:50*unite_base_u);
region_2       = signature_binaire2(51*unite_base_u:92*unite_base_u); % représente les chiffres de rang 8 à 13
garde_droite   = signature_binaire2(93*unite_base_u:95*unite_base_u);

% Construction de signature théorique S_th de longueur 7
% 3 familles de 10 configurations
element_A = [
                    [1 1 1 0 0 1 0]; % 0
                    [1 1 0 0 1 1 0]; % 1
                    [1 1 0 1 1 0 0]; % 2
                    [1 0 0 0 0 1 0]; % 3
                    [1 0 1 1 1 0 0]; % 4
                    [1 0 0 1 1 1 0]; % 5
                    [1 0 1 0 0 0 0]; % 6
                    [1 0 0 0 1 0 0]; % 7
                    [1 0 0 1 0 0 0]; % 8
                    [1 1 1 0 1 0 0]; % 9
                    ];
element_B = [
                    [1 0 1 1 0 0 0]; % 0
                    [1 0 0 1 1 0 0]; % 1
                    [1 1 0 0 1 0 0]; % 2
                    [1 0 1 1 1 1 0]; % 3
                    [1 1 0 0 0 1 0]; % 4
                    [1 0 0 0 1 1 0]; % 5
                    [1 1 1 1 0 1 0]; % 6
                    [1 1 0 1 1 1 0]; % 7
                    [1 1 1 0 1 1 0]; % 8
                    [1 1 0 1 0 0 0]; % 9
                    ];
element_C = [
                    [0 0 0 1 1 0 1]; % 0
                    [0 0 1 1 0 0 1]; % 1
                    [0 0 1 0 0 1 1]; % 2
                    [0 1 1 1 1 0 1]; % 3
                    [0 1 0 0 0 1 1]; % 4
                    [0 1 1 0 0 0 1]; % 5
                    [0 1 0 1 1 1 1]; % 6
                    [0 1 1 1 0 1 1]; % 7
                    [0 1 1 0 1 1 1]; % 8
                    [0 0 0 1 0 1 1]; % 9
                    ];
                
% Dilatation de la signature théorique en fonction de l’unité de base 
n_ligne = 10; 
n_colone = 7;
dilat_element_A = dilater_signature_th(element_A, n_ligne, n_colone, unite_base_u);
dilat_element_B = dilater_signature_th(element_B, n_ligne, n_colone, unite_base_u);
dilat_element_C = dilater_signature_th(element_C, n_ligne, n_colone, unite_base_u);

signature_partielle = double([region_1; region_2]);
mr_A = zeros(1,n_ligne);
mr_B = zeros(1,n_ligne);
mr_C = zeros(1,n_ligne);
for i=1:n_ligne
    mr_A(i) = norm(dilat_element_A(i,:) - signature_partielle((i-1)*unite_base_u+1:i*unite_base_u), 1);
    mr_B(i) = norm(dilat_element_B(i,:) - signature_partielle((i-1)*unite_base_u+1:i*unite_base_u), 1);
    mr_C(i) = norm(dilat_element_C(i,:) - signature_partielle((i-1)*unite_base_u+1:i*unite_base_u), 1);
end

figure, 
subplot(131), plot(mr_A), title('mrA')
ylim([0 8])
subplot(132), plot(mr_C), title('mrB')
ylim([0 8])
subplot(133), plot(mr_B), title('mrC')
ylim([0 8])
% affichage des mesure de ressemblance de la signature théorique avec 
% la signature partielle observée par comptage du nombre de différences 
% mesure_ressemblance = norm();


% Identifier les chiffres du code-barres
% chiffres_CB = zeros(1,13);

%% Fonctions
function signature_out = dilater_signature_th(signature_in, n_ligne, n_colone, unite_base)
    dilater_1 = ones(1,unite_base);
    signature_out = zeros(10, 7*unite_base);
    for i=1:n_ligne
        for j=1:n_colone
            if signature_in(i,j) == 1
                signature_out(i, ((j-1)*unite_base)+1:j*unite_base) = dilater_1;
            end
        end
    end
end
