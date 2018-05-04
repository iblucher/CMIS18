close all; 
clear all; 
clc;

k = [10, 25, 50, 100, 200, 500, 1000];
dx = zeros(size(k,2), 1);
time = zeros(size(k, 2), 1);

for i = 1:size(k,2)
    tic;
    [x, e] = advection(k(i));
    dx(i) = x;
    time(i) = toc;
end

plot(dx, time);
xlabel('$\Delta t$', 'Interpreter', 'latex');
ylabel('Runtime in seconds');