% --- Linear elasticity FEM Theory summary --------------------------------
%
% The continous displacement field is given by interpolation of the discrete
% nodal displacement field
%
%  u(x) = N * Ue;                         % This is a 2-by-1
%
% where 
%
%   N = [ w_1 *I, w_2 *I, w_3 *I ];      % This is 2-by-6
%
% and
%
%  Ue = [Ux(T(k,:));  Uy(T(k,:))];        % This is 6-by-one 
%
% and linear interpolation is used (bary-centric coordinates) 
%
%   w_1  =  area(x,2,3)  / A(k)
%   w_2  =  area(3,1,x)  / A(k)
%   w_3  =  area(x,1,2)  / A(k)
%
% The Cauchy strain tensor is given by
%
%  epsilon(x) = S * u(x) = S * N * Ue  = B * Ue; % This is 3-by-1  
%
% where 
%
%   B = S * N;                           % This is 3-by-6
%
% and the differential operator is defined as
%
%     [ dx  0;
% S =   0  dy; 
%       dy dx ];                         % This is 3-by-2
%
% The Cauchy Stress tensor can by found by the equation of state (EOS)
% which defines a relation between stress and strain. For a neo-hookean
% material model this yields the relation
%
%  sigma =  D epsilon;                             % This is 3-by-1
%
% where
%
%  D =  E/(1-nu^2) *[  1 nu 0;
%                     nu  1 0;
%                      0  0 (1-nu)/2];  % This is 3-by-3
%
% Now by using the virtual work one ends up with a linear system for the
% unknown nodal displacements
%
%  Ke Ue = Fe 
%
% where
%
%  Ke = B' * D * B * Ae;         % This is 6-by-6
%
% To assemble the element-wise equations into one large simultaneous system
% one applies Newton's third law at the nodes.
%
%--------------------------------------------------------------------------

clear all;
close all;

% --- Create the mesh geometry that we will use, a retangular bar --------- 
% --- If one has DistMesh installed then one can out-comment the code
% --- below in order to experiment with the mesh generation
%
% fd= @(p) drectangle(p,-3,3,-1,1)
% [p, T] = distmesh2d(fd,@huniform,0.5,[-3,-3;3,3],[-3,-1; -3,1; 3,-1; 3,1]);
%
% X      = p ( : , 1 ) ;
% Y      = p ( : , 2 ) ;
% cntV   = size(X,1);
% V      = (1:cntV)';
% cntT   = size(T,1);
%
% 
% dx1 = X( T(:,3) ) - X( T(:,2) );
% dy1 = Y( T(:,3) ) - Y( T(:,2) );
% dx2 = X( T(:,1) ) - X( T(:,2) );
% dy2 = Y( T(:,1) ) - Y( T(:,2) );
% A =  (dx1.*dy2 - dx2.*dy1 ) ./ 2; % Compute triangle areas
%
% clear p dx1 dx2 dy1 dy2 fd;
%
% save('dist.mat');
%
% --- If one does not have DistMesh then one can simply load a precomputed
% --- mesh that we have prepared  
load('dist.mat');

Ke = zeros(6,cntT*6);        % Allocate array of element stiffness matrices
K  = zeros(cntV*2,cntV*2);   % Allocate global stiffness matrix
B = zeros(3, 6);

% Now compute element stiffness matrices
for e=1:cntT
  
  % Get triangle indices
  i = T(e,1);  j = T(e,2); k = T(e,3);
  
  % Get triangle coordinates
  xi = X(i); xj = X(j); xk = X(k);
  yi = Y(i); yj = Y(j); yk = Y(k);
  
   
  % Compute the D (elasticity) matrix
  %
  % http://en.wikipedia.org/wiki/Young%27s_modulus
  % http://en.wikipedia.org/wiki/Poisson%27s_ratio
  %
  % This is steel like parameters
  %
  E  = 69e9;    % Young modulus
  nu = 0.3;    % Poisson ration
  D  = E/(1-nu.^2) *[  1 nu 0; nu  1 0; 0  0 (1-nu)/2];
 
  % Do cross product to get the area
  V1 = [xj - xi, yj - yi];
  V2 = [xk - xi, yk - yi];
  Ae = 0.5 * (V1(1) * V2(2) - V1(2) * V2(1));
  
  delNx = [-(yk - yj)/(2*Ae), -(yi - yk)/(2*Ae), -(yj - yi)/(2*Ae)];
  delNy = [(xk - xj)/(2*Ae), (xi - xk)/(2*Ae), (xj - xi)/(2*Ae)];
  
  B(1, [1, 3, 5]) = delNx;
  B(2, [2, 4, 6]) = delNy;
  B(3, [1, 3, 5]) = delNy;
  B(3, [2, 4, 6]) = delNx;
  
  Ke( 1:6, (6*(e-1)+1) :(6*e)) =  B' * D * B * Ae;
  
end

% Assembly of  element matrix fe
fy = 0;
fx = 10e9;
ff = [fx; fy];
indices = find(X > 2.9);
fe = zeros(4, length(indices) - 1);

sortedY=Y(indices);
sortedY=sort(sortedY);
for v = 1:(length(indices) - 1)
    xi = X(indices(v));
    xj = X(indices(v + 1));
    yi= sortedY(v,1); 
    yj = sortedY(v+1,1);
    I = [1 0 1 0; 0 1 0 1];
    Le = sqrt((xj - xi)^2 + (yj - yi)^2);
    fe(:, v) = 0.5 * Le * I' * ff;
end
   
% Now do assembly process of global stiffness matrix
for e=1:cntT
  
  % Get global triangle vertex indices
  i = T(e,1);  j = T(e,2); k = T(e,3);
  
  % Column offset of the local element stiffness matrix into the Ke array
  coffset = + 6*(e-1);  
  
  % Local order of vertex coordinates is i_x, i_y, j_x j_y, k_x, and  k_y. 
  % This is how local vertex indices (1,2,3,..,6) are mapped to global vertex
  % indices
  gidx = [ i, cntV + i, j,  cntV + j,  k, cntV + k];

  % Now we can add the element stiffness matrix to the global matrix using
  % our local-global vertex index mapping.
  K(gidx, gidx) = K(gidx, gidx) + Ke(1:6, coffset+(1:6)); 
  
end

figure(1);
clf;  
hold on;  
spy(K);
title('Fill pattern of stiffness matrix');
ylabel('Row index');
xlabel('Column index');
hold off;

figure(2);
clf;  
hold on;  
plot(sort( eig(K) ), 'r-', 'LineWidth', 2);
title('Eigenvalues of stiffness matrix');
xlabel('Eigenvalue Index')
ylabel('Value');
hold off;
axis tight;

% Create an  external nodal force vector with a prescribed load at the
% right hand side of the bar

load = -10e10;
f = zeros(2*cntV,1);
indices = find(X>2.9);
%f( indices + cntV ) =  load;
% Assembly process for f matrix
for a = 1:(length(fe) - 1)
    i = indices(a);
    j = indices(a + 1);
    f(i) = f(i) + fe(1, a);
    f(i + cntV) = f(i + cntV) + fe(2, a);
    f(j) = f(j) + fe(3, a);
    f(j + cntV) = f(j + cntV) + fe(4, a);
end

% Apply Direchlet boundary conditions to displacement field for the left
% vertices.
indices = find(X<-2.9);
indices = [indices; indices+cntV];
values  = zeros(size(indices));

F = setdiff( 1:cntV*2, indices); 

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

% Compute displacement
u          = zeros(size(2*cntV,1));
u(indices) = values;
u(F)       = K(F,F) \ f(F);

% Compute deformed coordinates
x = X + u(1:cntV)';
y = Y + u(cntV+1:end)';
indices = find(X>2.9);
sortedY=y(indices);
sortedY=sort(sortedY);
corner_y = sortedY(1); %CORNER NODE

figure(3);
clf;
hold on;
triplot(T,X,Y,'b');
triplot(T,x,y,'r');
title('The computational mesh');
xlabel('x');
ylabel('y');
axis equal;
hold off;