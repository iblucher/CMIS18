function [dt, error] = advection(k)

% Initializing variables
n = 32;        % length of grid
%k = 100;        % number of steps to reach max time
L = 10;         % used for calculating grid interval
T = 5;          % maximum time 
dt = T/k;       % time step
dx = L/n;       % grid step size
dy = dx;

% Initialize peaks as scalar field
[X, Y] = meshgrid(-L/2:dx:L/2, -L/2:dy:L/2);
Z = peaks(X, Y);
figure(1);
surf(Z);

% Backtrack coordinates
nX = X;
nY = Y;

% Update equations
uX = Y * dt;
uY = (-X) * dt;
nZ = Z;

% Energy of scalar field
E0 = sum(sum(dx * dy * nZ));

% Time dependent loop
for t = 1:k
    t
    nX = X - uX;
    nY = Y - uY;
    
    %interpolate
    nZ = interp2(X, Y, nZ, nX(:), nY(:), 'linear', 0);
    nZ = reshape(nZ, [n + 1, n + 1]);
    figure(1);
    surf(X, Y, nZ);
    zlim([-10, 10]);
    %pause(0.05)
end

nE = sum(sum(dx * dy * nZ));
error = abs(nE - E0);
end