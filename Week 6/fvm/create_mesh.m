%--- Create the mesh geometry that we will use, a retangular bar --------- 
%--- If one has DistMesh installed then one can out-comment the code
%--- below in order to experiment with the mesh generation

fd= @(p) drectangle(p,-3,3,-1,1)
[p, T] = distmesh2d(fd,@huniform,0.15,[-3,-3;3,3],[-3,-1; -3,1; 3,-1; 3,1]);

X      = p ( : , 1 ) ;
Y      = p ( : , 2 ) ;
cntV   = size(X,1);
V      = (1:cntV)';
cntT   = size(T,1);


dx1 = X( T(:,3) ) - X( T(:,2) );
dy1 = Y( T(:,3) ) - Y( T(:,2) );
dx2 = X( T(:,1) ) - X( T(:,2) );
dy2 = Y( T(:,1) ) - Y( T(:,2) );
A =  (dx1.*dy2 - dx2.*dy1 ) ./ 2; % Compute triangle areas

clear p dx1 dx2 dy1 dy2 fd;

save('dist.mat');