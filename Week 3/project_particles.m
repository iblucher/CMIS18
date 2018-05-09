function [X, Y] = project_particles(phi, X, Y)
    [M, N] = size(phi);
    [GX, GY] = meshgrid(1:N, 1:M);
    d = interp2(GX, GY, phi, X, Y);
    [dx, dy] = gradient(phi);
    nx = interp2(GX, GY, dx, X, Y);
    ny = interp2(GX, GY, dy, X, Y);
    dx = d .* nx;
    dy = d .* ny;
    X(d>0) = X(d>0) - dx(d>0);
    Y(d>0) = Y(d>0) - dy(d>0);
end