clear all;
close all;

dx = 0.1;
x  = 1:dx:2;
N  = length(x);
N
Ke = zeros(2,(N-1)*2);  
K  = zeros(N,N);
f  = zeros(N,1);

for e=1:(N-1)  
  i = e;
  j = i+1;
 
  Ke( 1:2, (2*(e-1)+1) :(2*e)) =  1/dx * [1, -1; -1, 1];
end

for e=1:(N-1)  
  i = e;
  j = i+1;
  offset = 2*(e-1);
  
  K(i, i) = K(i, i) + Ke(1, 1 + offset);
  K(i, j) = K(i, j) + Ke(1, 2 + offset);
  K(j, i) = K(j, i) + Ke(2, 1 + offset);
  K(j, j) = K(j, j) + Ke(2, 2 + offset);
  
end

figure(1);
clf;  
hold on;  
spy(K);
title('Fill pattern of matrix');
ylabel('Row index');
xlabel('Column index');
hold off;

figure(2);
clf;  
hold on;  
plot(sort( eig(K) ), 'r-', 'LineWidth', 2);
title('Eigenvalues of matrix');
xlabel('Eigenvalue Index')
ylabel('Value');
hold off;
axis tight;

% Apply boundary conditions
a = 1;
b = 2;
indices = [1; N];
values  = [a; b];
F       = setdiff( 1:N, indices);

for i = 1: length(indices)
  index = indices(i);
  value = values(i);  
  
  f(index) = value;
  for j = indices(1) + 1:indices(2) - 1
      f(j) = f(j) - K(j, index) * value;
  end
  
  K(index, :) = 0;
  K(:, index) = 0;
  K(index, index) = 1;

end

y = zeros(size(f));
y(indices) = values;
y(F)       = K(F,F) \ f(F);

figure(3);
clf;
hold on;
plot(x,y,'b');
title('The approximate solution');
xlabel('x');
ylabel('y');
axis equal;
hold off;
