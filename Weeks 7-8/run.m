clear all;
close all;
clc;

% The triangle mesh
load('data.mat');  % Load mesh

[TR CVs Bmask] = create_control_volumes(T, X, Y);

figure(1);
clf;
draw_control_volumes(TR, CVs, 0.05);
axis equal;
title('Computational mesh')

% t = % Surface traction field
% 
% rho = % Material density field
% 
% b = % Body force density
% 
% [lam, mu] = % The Lamé parameters (elasticity)
% 
% dt = % The time step
% 
% simulator(CVs, T, X, Y, cntT, t, rho, b, lam, mu, dt);
simulator(CVs, T, X, Y, cntT);

