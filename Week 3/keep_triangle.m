function keep_triangle(n, T, i, phi)

[M, N] = size(phi);
[GX, GY] = meshgrid(1:N, 1:M);

while count < n:
    % generate random point
    [X, Y] = barycentricToCartesian(T, i, rand(1, 3));
    
    % check if in polygon
    
    % interpolation
end

end