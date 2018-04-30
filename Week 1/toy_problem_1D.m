% 1D Toy Problem for HW1 (CMIS18)

xa = 0;
xb = 1;

n = 4; % number of points sampled on interval [0, 1]
dx = 1/(n+1);
N = n + 2; % number of actual domain nodes (includes 0 and 1)
N_g = N + 2; % number of nodes including ghost nodes 

% Create mesh interval
mesh_interval = (xa:dx:xb);
mesh_interval

% Draw computational mesh
figure(1);
plot(mesh_interval, 1, '-.ro');

% Assemble matrix system
% Compute coefficients of domain nodes
c0 = -2.0/dx^2;
c1 = 1.0/dx^2;

A = zeros(N_g);
for i = 2:N_g-1
    k = i;
    A(k, k) = c0;
    A(k, k - 1) = c1;
    A(k, k + 1) = c1;
end

% Compute coefficients of ghost nodes
A(1, 1) = 1;
A(1, 2) = -1;
A(8, 8) = 1;
A(8, 7) = -1;

% Plot matrix structure
f_type = 'Times';
f_size = 12;

figure(2);
imagesc(A);
grid on;
colorbar;
title('Assembly matrix structure', 'FontSize', f_size, 'FontName',f_type );
xlabel('Columns', 'FontSize', f_size, 'FontName', f_type); 
ylabel('Rows', 'FontSize', f_size, 'FontName', f_type);

% Compute matrix eigenvalues
e = eig(A);
