% Initialize variables
dx = dy = 1;
kappa = 2;
n = 4;

% Compute f
f = zeros(n, n);
for i = 1:n
    for j = 1:n
        f(i, j) = i + j;
    end
end

% Compute coefficients
c0 = -kappa^2 - 2/dx^2 - 2/dy^2;
c1 = 1/dx^2;
c2 = 1/dy^2;

% Initialize domain coefficient matrix
N = (n + 2)^2;
A = zeros(N);
for j = 2:N-1
    for i = 2:N-1
        k = j*n + i;
        A(k, k) = c0;
        A(k, k - 1) = c1;
        A(k, k + 1) = c1;
        A(k, k - n) = c2;
        A(k, k + n) = c2;
    end
end

% 
    

% compute A
           
            
