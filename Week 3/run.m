clear all;
close all;

I = read_bw( 'EG_WEB_logo.jpg');  % Read black and white image from file
phi = bw2phi(I);

imagesc(phi);



