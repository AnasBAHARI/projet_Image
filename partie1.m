clear all ;
close all ;
clc;

%% Lecture

img = imread('code_barre2.png');
[h, w, c] = size(img);
fprintf("L'image est une matrice de %s entre %d et %d\n", class(img(10,10,2)), min(img(:)), max(img(:)));
fprintf("Elle est de dimension: %d*%d*%d\n", h, w, c);
figure, imshow(img);

%[x, y] = ginput(2);
%fprintf('les points choisies sont (%i, %i) et (%i, %i)\n', x(1), y(1), x(2), y(2));
x = [4.409999999999999e+02;79];
y = [56.999999999999940;3.869999999999999e+02];
Ln = sqrt( (x(1)-x(2))^2 + (y(1)-y(2))^2 ); % longueur du segment

U = 2*ceil(Ln); % Nombre de point a prendre

Mu = zeros(U, 2);
Mu(1,:) = [x(1) y(1)];
Mu(U-1,:) = [x(2) y(2)];
for u=2:U-2
    Mu(u,:) = Mu(1,:) + (u/(U-1))*(Mu(U-1,:) - Mu(1,:));
end

%hold on
%plot(Mu(:,1), Mu(:,2), '*');

for u=1:U-1
    pt_img = sum(img( ceil(Mu(u,1)), ceil(Mu(u,2)), :));
    vect_afficher(u) = pt_img;
end
figure,
plot(vect_afficher(1,:))

