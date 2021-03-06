clear all;
close all;

load('data.mat');  % Load mesh

Ke = zeros(3,cntT*3);
K  = zeros(cntV,cntV); 
f  = zeros(cntV,1);

for e=1:cntT
  
  % Get triangle indices
  i = T(e,1);  j = T(e,2); k = T(e,3);
  
  % Get triangle coordinates
  xi = X(i); xj = X(j); xk = X(k);
  yi = Y(i); yj = Y(j); yk = Y(k);
  
  % Do cross product to get the area
  V1 = [xj - xi, yj - yi];
  V2 = [xk - xi, yk - yi];
  Ae = 0.5 * (V1(1) * V2(2) - V1(2) * V2(1));
  
  % Compute partial derivatives and assemble Ke
  delNx = [-(yk - yj)/(2*Ae), -(yi - yk)/(2*Ae), -(yj - yi)/(2*Ae)];
  delNy = [(xk - xj)/(2*Ae), (xi - xk)/(2*Ae), (xj - xi)/(2*Ae)];
  Ke( 1:3, (3*(e-1)+1) :(3*e)) =  (delNx' * delNx + delNy' * delNy) * Ae;  
end

for e=1:cntT
  
  % Get global triangle vertex indices
  i = T(e,1);  j = T(e,2); k = T(e,3);
  
  % Column offset of the local element matrix into the Ke array
  coffset = + 3*(e-1);  
  
  % Local order of vertex coordinates is i j and k. 
  % This is how local vertex indices (1,2,3) are mapped to global vertex
  % indices
  gidx = [ i, j,  k];

  for ii = 1:3
      for jj = 1:3
          K(gidx(ii), gidx(jj)) = K(gidx(ii), gidx(jj)) + Ke(ii, jj + coffset);
      end
  end
  
end

% figure(1);
% clf;  
% hold on;  
% spy(K);
% title('Fill pattern of stiffness matrix');
% ylabel('Row index');
% xlabel('Column index');
% hold off;

% figure(2);
% clf;  
% hold on;  
% plot(sort( eig(K) ), 'r-', 'LineWidth', 2);
% title('Eigenvalues of stiffness matrix');
% xlabel('Eigenvalue Index')
% ylabel('Value');
% hold off;
% axis tight;

% Apply boundary conditions
left_edge   = find(X<-2.9);
right_edge  = find(X>2.9);
middle_edge = find(X<.1 & X>-0.1);
a = 1;
b = 2;
%c = 3;
%indices    = [left_edge; middle_edge; right_edge];
indices    = [left_edge; right_edge];
%values     = [ones(size(left_edge))*a ; ones(size(middle_edge))*c ;  ones(size(right_edge))*b];
values     = [ones(size(left_edge))*a;  ones(size(right_edge))*b];


u = zeros(size(X));
F = setdiff( 1:cntV, indices);

for i = 1: length(indices)
  index = indices(i);
  value = values(i);

  f(index) = value;
  for j = 1:cntV
      if ismember(j, indices)
          continue;
      end
      f(j) = f(j) - K(j, index) * value;
  end
  
  K(index, :) = 0;
  K(:, index) = 0;
  K(index, index) = 1;
  
end

% Compute solution
u(indices) = values;
u(F)       = K(F,F) \ f(F);

% Visualize solution
figure(3);
clf;
pX = [X(T(:,1)) X(T(:,2)) X(T(:,3)) ]';
pY = [Y(T(:,1)) Y(T(:,2)) Y(T(:,3)) ]';
pZ = [u(T(:,1)) u(T(:,2)) u(T(:,3)) ]';
patch( pX, pY, pZ, pZ );
hold on;
view(3)
colorbar
shading interp;
title('The Solution');
xlabel('x');
ylabel('y');
zlabel('u');
axis equal;
trisurf(T, X, Y, u);
hold off;