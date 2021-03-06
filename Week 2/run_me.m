close all;
clear all;
clc;

U = double(im2bw(imread('example.bmp')));
phi = bw2phi(U);

fig1 = figure(1);
imagesc(phi); colormap(gray); axis image; colorbar
hold on; title('Input');contour(phi,[0,0],'r'); hold off
drawnow;

dx = [0.5, 0.75, 1, 1.25, 1.5, 2];
time = zeros(size(dx, 2), 1);
for i = 1:size(dx, 2)
    tic;
    psi = mean_curvature_flow(phi, 100.0, dx(i)); 
    time(i) = toc;
end

fig2 = figure(3);
imagesc(psi); colormap(gray); axis image; colorbar
hold on; title('Output');contour(phi,[0,0],'r'); hold off
drawnow;


