% Initialize variables
dx = 1;
dy = 1;
kappa = 2;
n = 4;

% Compute f
f = zeros(n, n);
for i = 1:n+2
    for j = 1:n+2
        f(i, j) = i + j;
    end
end

imagesc(f)

% Compute A
% Compute coefficients
c0 = -kappa^2 - 2/dx^2 - 2/dy^2;
c1 = 1/dx^2;
c2 = 1/dy^2;

% Initialize domain nodes coefficient
N = (n + 2)^2;
A = zeros(N);
for j = 2:n+1
    for i = 2:n+1
        k = j*n + i;
        A(k, k) = c0;
        A(k, k - 1) = c1;
        A(k, k + 1) = c1;
        A(k, k - n) = c2;
        A(k, k + n) = c2;
    end
end

% Initialize ghost nodes coefficients
for i = 2:n+1
    k1 = 1*n + i;
    k2 = 6*n + i;
    A(k1, k1) = 1;
    A(k1, k1 + n) = -1;
    A(k2, k2) = 1;
    A(k2, k2 - n) = -1;
end

for j = 2:n+1
    k3 = j*n + 1;
    k4 = j*n + 6;
    A(k3, k3) = 1;
    A(k3, k3 + 1) = -1;
    A(k4, k4) = 1;
    A(k4, k4 - 1) = -1;
end

f_new = reshape(f, [(n+2)^2, 1]);

% Compute solution for Au = f

u = A \ f_new;


           
            
