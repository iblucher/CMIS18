clear all;
close all;

I = read_bw( 'EG_WEB_logo.jpg');  % Read black and white image from file
phi = bw2phi(I);

figure(1);
imagesc(phi);
colormap bone;
[M, N] =size(phi);
K=1000;
Y = rand(1,K)*(M-6)+3;
X = rand(1,K)*(N-6)+3;
figure(2);
plot(X,Y,"r+")
xlabel('x coordinates');
ylabel('y coordinates');
%axis([0 120 0 120]);
[GX,GY] = meshgrid(1:N,1:M);

figure(3);

for i =1:5
%[X,Y]=spring(X,Y);
[X,Y]=project(X,Y,phi,GX,GY);
plot(X,Y,"r+")
xlabel('x coordinates');
ylabel('y coordinates');
%pause(1);
end
%%
T = delaunay(X, Y);
figure(4);
triplot(T, X, Y);
xlabel('x coordinates');
ylabel('y coordinates');

legitTris=delMesh(X,Y,phi);
%%
triplot(legitTris,X,Y);
xlabel('x coordinates');
ylabel('y coordinates');
patch('Faces',legitTris,'Vertices',[X;Y]','FaceColor','red')