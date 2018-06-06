function simulator(CVs, T, X, Y, cntT)
    [N, ~] = size(CVs);
    t = 0;

    % Computing the original deformation gradient D0 
    % These are material coordinates, so they are going to be constant to
    % time changes
    D0 = zeros(2, cntT * 2);
    X0 = X;
    Y0 = Y;
    for tri = 1:cntT
       i = T(tri, 1);
       j = T(tri, 2);
       k = T(tri, 3);
       
       offset = (tri - 1) * 2 + 1;
       Gjix = abs(X0(j) - X0(i));
       Gjiy = abs(Y0(j) - Y0(i));
       Gkix = abs(X0(k) - X0(i));
       Gkiy = abs(Y0(k) - Y0(i));
       
       temp = [Gjix Gkix; Gjiy Gkiy];
       D0e = pinv(temp);
       D0(:, offset:offset + 1) = D0e;
    end
    
    
    D = zeros(2, cntT * 2);
    F = zeros(2, cntT * 2);
    E = zeros(2, cntT * 2);
    S = zeros(2, cntT * 2);
    P = zeros(2, cntT * 2);
    while t <= 1
        % Iterating on the triangle mesh
        for tri = 1:cntT
            % compute D
            i = T(tri, 1);
            j = T(tri, 2);
            k = T(tri, 3);
       
            offset = (tri - 1) * 2 + 1;
            Gjix = abs(X(j) - X(i));
            Gjiy = abs(Y(j) - Y(i));
            Gkix = abs(X(k) - X(i));
            Gkiy = abs(Y(k) - Y(i));
       
            temp = [Gjix Gkix; Gjiy Gkiy];
            D(:, offset:offset + 1) = temp;
            
            % compute F
            % compute deformation gradients Fe   
            F(:, offset:offset + 1) = D(:, offset:offset + 1) * D0(:, offset:offset + 1);
            
            % compute E
            % compute green strain tensors Ee
            Fe = F(:, offset:offset + 1); 
            E(:, offset:offset + 1) = 0.5 * (Fe' * Fe - eye(2));
            
            % compute S
            % compute 2nd piola-kirchhoff stress tensor Se
            Ee = E(:, offset:offset + 1);
            S(:, offset:offset + 1) = lambda * trace(Ee) * eye(2) + 2 * mu * Ee;
            
            % compute P
            % compute 1st piola-kirchhoff stress tensor Pe
            Se = S(:, offset:offset + 1);
            P(:, offset:offset + 1) = Fe * Se;
        end
 
        
        for i = 1:N
            CV = CVs{i};
            
             % sum all fe for this i
             % compute ftotal for this i
             
             % update v for this i
             % update x for this i
        end
        t = t + 1;
    end
    
end