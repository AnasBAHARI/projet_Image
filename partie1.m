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
hold on, plot(k_max*ones(nbr_pt));
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
hold on, plot(k_max*ones(nbr_pt2));
title("Signature de l'image entre les points utiles avec le seuil");

signature_binaire2 = signature_img2 > k_max;
figure, plot(signature_binaire2);
ylim([-0.25 1.25]);
title("Signature Binaire de l'image entre les points utiles");


