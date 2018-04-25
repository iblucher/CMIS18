% Initialize variables
dx = 1;
dy = 1;
kappa = 2;
n = 4;
N = n + 2;

% Compute f
f = zeros(N);
for i = 2:N-1
    for j = 2:N-1
        f(i, j) = i + j;
    end
end

imagesc(f)

% Compute A
% Compute coefficients
c0 = -kappa^2 - 2.0/dx^2 - 2.0/dy^2;
c1 = 1.0/dx^2;
c2 = 1.0/dy^2;

% Initialize domain nodes coefficients
A = zeros(N^2);
for j = 2:N-1
    for i = 2:N-1
        k = (j-1)*N + i;
        A(k, k) = c0;
        A(k, k - 1) = c1;
        A(k, k + 1) = c1;
        A(k, k - N) = c2;
        A(k, k + N) = c2;
    end
end

% Initialize ghost nodes coefficients
for i = 2:N-1
    j1 = 1;
    j2 = 6;
    k1 = (j1 - 1) * N + i;
    k2 = (j2 - 1) * N + i;
    A(k1, k1) = 1;
    A(k1, k1 + N) = -1;
    A(k2, k2) = 1;
    A(k2, k2 - N) = -1;
end

for j = 2:N-1
    i1 = 1;
    i2 = 6;
    k3 = (j-1)*N + i1;
    k4 = (j-1)*N + i2;
    A(k3, k3) = 1;
    A(k3, k3 + 1) = -1;
    A(k4, k4) = 1;
    A(k4, k4 - 1) = -1;
end

A(1, 1) = 10;
A(N, N) = 10;
A(N*N, N*N) = 10;
A(N*(N-1) + 1, N*(N-1) + 1) = 10;

f_new = reshape(f, [(N)^2, 1]);

imagesc(A)

% Compute solution for Au = f

u = A \ f_new;


% Compute eigenvalues
e = eig(A);


           
            
