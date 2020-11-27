clear all ;
close all ;
clc;

%% Lecture

img = imread('code_barre.png');
[h, w, c] = size(img);
fprintf("L'image est une matrice de %s entre %d et %d\n", class(img(10,10,2)), min(img(:)), max(img(:)));
fprintf("Elle est de dimension: %d*%d*%d\n", h, w, c);
%figure, imshow(img);

% [x, y] = ginput(2);
% fprintf('les points choisies sont (%i, %i) et (%i, %i)\n', x(1), y(1), x(2), y(2));
x = [885;531];
y = [175;5050000000000001e+02];
Ln = norm([x(1)-x(2), y(1)-y(2)], 2);

U = 2*ceil(Ln) % Nombre de point a prendre

