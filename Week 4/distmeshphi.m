I = read_bw( 'EG_WEB_logo.jpg');  % Read black and white image from file
phi = bw2phi(I);
[M, N] =size(phi);
[GX,GY] = meshgrid(1:N,1:M);
fd2=@(p) interp2(GX,GY,phi,p(:,1),p(:,2 ));
[p,t]=distmesh2d(fd2,@huniform,3,[1,1;170,110],[]);