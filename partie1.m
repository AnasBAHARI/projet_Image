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
% x = [4.489999999999999e+02;77.999999999999990];
% y = [43.999999999999940;4.009999999999999e+02];

%% Lancer Aléatoire d'un rayon
Ln = sqrt( (x(1)-x(2))^2 + (y(1)-y(2))^2 ); % longueur du segment

nbr_pt = 2*ceil(Ln); % Nombre de point a prendre on respectant le Crtr2Shannon

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

%% Extraction de la signature le long d’un rayon
figure, hist(signature_img, 256);

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
hold on,
seuil = k_max*ones(length(signature_img));
plot(seuil);

signature_binaire = (signature_img>k_max);
figure, plot(signature_binaire);
ylim([-0.5 1.5])
signature_binaire = signature_binaire<1;
figure, plot(signature_binaire);
ylim([-0.5 1.5])





