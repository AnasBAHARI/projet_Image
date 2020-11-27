clear all;
close all;
clc;

img = imread("barcode1D.jpg");
img2 = imread("rotated1DBarcode.jpg");


figure, imshow(img);
figure, imshow(img2);