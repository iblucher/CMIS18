close all;
clear all;

U = double(im2bw(imread('example.bmp')));
phi = bw2phi(U);

fig1 = figure(1);
imagesc(phi); colormap(gray); axis image; colorbar
hold on; title('Input');contour(phi,[0,0],'r'); hold off
drawnow;

phi = mean_curvature_flow(phi, 1000.0); 

fig2 = figure(3);
imagesc(phi); colormap(gray); axis image; colorbar
hold on; title('Output');contour(phi,[0,0],'r'); hold off
drawnow;


