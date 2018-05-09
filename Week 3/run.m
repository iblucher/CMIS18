clear all;
close all;

I = read_bw( 'EG_WEB_logo.jpg');  % Read black and white image from file
phi = bw2phi(I);

figure(1);
imagesc(phi);

figure(2);
[X, Y] = random_particles(phi, 1000);
plot(X, Y, 'r+');

figure(3);
[nX, nY] = project_particles(phi, X, Y);
plot(nX, nY, 'r+');
hold on;

T = delaunay(nX, nY);
triplot(T, nX, nY);
hold off;
