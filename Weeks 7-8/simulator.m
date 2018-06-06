function simulator(CVs, T, X, Y, cntT)
    [N, ~] = size(CVs);
    t = 0;

    % Computing the original deformation gradient D0 
    % These are material coordinates, so they are going to be constant to
    % time changes
    Gji = zeros(cntT, 2);
    Gki = zeros(cntT, 2);
    for tri = 1:cntT
       Gji(tri, 1) = abs(X(T(tri, 2)) - X(T(tri, 1)));
       Gji(tri, 2) = abs(Y(T(tri, 2)) - Y(T(tri, 1)));
       Gki(tri, 1) = abs(X(T(tri, 3)) - X(T(tri, 1)));
       Gki(tri, 2) = abs(Y(T(tri, 3)) - Y(T(tri, 1)));
    end
    
     
    while t <= 1
        % Iterating on the triangle mesh
        for tri = 1:cntT
            % compute deformation gradients Fe    
            % compute green strain tensors Ee
            % compute 2nd piola-kirchhoff stress tensor Se
            % compute 1st piola-kirchhoff stress tensor Pe
            % compute elastic force fe 
        
        
        for i = 1:N
            CV = CVs{i};
            
             E = length(CV.I); % Get the number of edge on the boundary of this CV
             for e=1:E
                 % compute deformation gradients Fe
                 
                 % compute green strain tensors Ee
                 % compute 2nd piola-kirchhoff stress tensor Se
                 % compute 1st piola-kirchhoff stress tensor Pe
                 
                 % compute elastic force fe 
             end
             
             % sum all fe for this i
             % compute ftotal for this i
             
             % update v for this i
             % update x for this i
        end
        t = t + dt;
    end
end