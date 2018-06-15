clear all;
close all;
clc;

% The triangle mesh
load('dist.mat');  % Load mesh

[TR CVs Bmask] = create_control_volumes(T, X, Y);

%figure(1);
%clf;
%draw_control_volumes(TR, CVs, 0.05);
%axis equal;
%title('Computational mesh');

gravity = -9.82;

rho = 1050;
b = [0; gravity * rho];
%b = [0; 0];

traction1 = [0; -1e5];
traction2 = [0; -1e5];

% The Lamé parameters (elasticity)
E = 0.01e9;
nu = 0.48; 

mu = E / (2 * (1 + nu));
lambda = (E * nu) / ((1 + nu) * (1 - 2 * nu));

%lambda = 1000;
%mu = 0.31;

% Time step
dt = 4000/ E;

X0 = X;
Y0 = Y;

Diri = find(X0 > -0.05 & X0 < 0.05);

%X = 1.1 * (X + 3) - 3;

simulator(CVs, T, X, Y, X0, Y0, cntT, lambda, mu, rho, b, dt, traction1, traction2, Diri);

