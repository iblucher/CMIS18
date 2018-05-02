function new_phi = ghost_nodes(phi)

    l = size(phi, 1);
    new_phi = zeros(l + 2);

    % Domain values to the center of new phi
    new_phi(2:l + 1, 2:l + 1) = phi;
    
    % Ghost rows and columns
    new_phi(2:l + 1, 1) = phi(1:l, 1);
    new_phi(2:l + 1, l + 2) = phi(1:l, l);
    new_phi(1, 2:l + 1) = phi(1, 1:l);
    new_phi(l + 2, 2:l + 1) = phi(l, 1:l);
end