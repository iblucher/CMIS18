function [X, Y] =  random_particles(phi, K)
    [M, N] = size(phi);
    Y = rand(1, K) * (M - 6) + 3;
    X = rand(1, K) * (N - 6) + 3;
end