function simulator(CVs, T, X, Y, X0, Y0, cntT, lambda, mu, rho, b0, dt, traction, traction2, Diri)
    [N, ~] = size(CVs);
    t = 0;
    velocities = zeros(2, N);
    count = 0;
    
    % find boundary edges 
    TR2 = triangulation(T,X,Y);
    es = edges(TR2);
    es1 = es(:,1);
    es2 = es(:,2);
    bound = es(X(es1) > 2.99 & X(es2) > 2.99,:);
    bound2 = es(X(es1) < -2.99 & X(es2) < -2.99,:);
    
    
    % Computing the original deformation gradient D0 
    % These are material coordinates, so they are going to be constant to
    % time changes
    D0 = zeros(2, cntT * 2);
    %X0 = X;
    %Y0 = Y;
    for tri = 1:cntT
       i = T(tri, 1);
       j = T(tri, 2);
       k = T(tri, 3);
       
       offset = (tri - 1) * 2 + 1;
       Gjix = X0(j) - X0(i);
       Gjiy = Y0(j) - Y0(i);
       Gkix = X0(k) - X0(i);
       Gkiy = Y0(k) - Y0(i);
       
       temp = [Gjix Gkix; Gjiy Gkiy];
       D0e = pinv(temp);
       D0(:, offset:offset + 1) = D0e;
    end

    D = zeros(2, cntT * 2);
    F = zeros(2, cntT * 2);
    E = zeros(2, cntT * 2);
    S = zeros(2, cntT * 2);
    P = zeros(2, cntT * 2);
    while t <= 1000
        X(Diri) = X0(Diri);
        Y(Diri) = Y0(Diri);
        velocities(:, Diri) = 0;
        
        % Iterating on the triangle mesh
        for tri = 1:cntT
            % compute D
            i = T(tri, 1);
            j = T(tri, 2);
            k = T(tri, 3);
       
            offset = (tri - 1) * 2 + 1;
            Gjix = X(j) - X(i);
            Gjiy = Y(j) - Y(i);
            Gkix = X(k) - X(i);
            Gkiy = Y(k) - Y(i);
       
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
        
        f_elastic = zeros(2, N);
        f_ext = zeros(2, N);
        area_cvs = zeros(1, N);
        f_total = zeros(2, N);
        ft = zeros(2, N);
        for cv = 1:N
            CV = CVs{cv};
            fe = 0;
            Ae = 0;
            inds = CV.indices;
            
            for e = 1:length(inds)
                tri_index = inds(e);
                a = T(tri_index, 1);
                b = T(tri_index, 2);
                c = T(tri_index, 3);
                center = CV.I(e);
                [i, j, k] = find_vertex_order(center, a, b, c);
                
                % access Pe
                Pe = P(:, tri_index * 2 - 1:tri_index * 2);
                
                % compute lengths
                xji = X0(j) - X0(i);
                yji = Y0(j) - Y0(i);
                xki = X0(k) - X0(i);
                yki = Y0(k) - Y0(i);
                Lji = sqrt((xji)^ 2 + (yji)^2);
                Lki = sqrt((xki)^ 2 + (yki)^2);
                
                % compute normal vectors
                Vji = [xji; yji];
                Vki = [xki; yki];
                Nji = [Vji(2); -Vji(1)]; % vector orientations may be flipped
                Nki = [-Vki(2); Vki(1)];
                
                Nji = Nji/abs(Lji);
                Nki = Nki/abs(Lki);
       
                % sum all fe for this i
                fe = fe + (- 0.5 * Pe * Nji * Lji - 0.5 * Pe * Nki * Lki);
                %fe = [100; 100];
                
                % compute area of triangles and sum to get area of control
                % volume
                Ae = Ae + 1/3 * abs(0.5 * (Vji(1) * Vki(2) - Vji(2) * Vki(1)));
                
            end
            f_elastic(:, cv) = fe;
            f_ext(:, cv) = Ae * b0;
            area_cvs(cv) = Ae;
             
        end
        
        % computing traction forces
        for e = 1:(size(bound,1))
            xi = X(bound(e,1));
            xj = X(bound(e,2));
            yi = Y(bound(e,1));
            yj = Y(bound(e,2));
            fl = traction;
            l_half = sqrt((xi - xj)^2 + (yi - yj)^2);
            l_half = l_half/2;
            ft(:,bound(e,1)) = ft(:,bound(e,1)) + fl * l_half;
            ft(:,bound(e,2)) = ft(:,bound(e,2)) + fl * l_half;  
        end

        for e = 1:(size(bound2,1))
            xi = X(bound2(e,1));
            xj = X(bound2(e,2));
            yi = Y(bound2(e,1));
            yj = Y(bound2(e,2));
            fl = traction2;
            l_half = sqrt((xi - xj)^2 + (yi - yj)^2);
            l_half = l_half/2;
            ft(:,bound2(e,1)) = ft(:,bound2(e,1)) + fl * l_half;
            ft(:,bound2(e,2)) = ft(:,bound2(e,2)) + fl * l_half;  
        end
        
        for cv = 1:N
            f_total(:, cv) = f_ext(:, cv) + f_elastic(:, cv) + ft(:, cv);
        end
        
        for cv = 1:N
            CV = CVs{cv};
            index = CV.I;
            index = index(1);
            
            % velocity updates
            mi = area_cvs(cv) * rho;
            f = f_total(:, cv);
            velocities(:, cv) = velocities(:, cv) + dt/mi * f;
            v = velocities(:, cv);
            
            % position updates
            X(index) = X(index) + dt * v(1);
            Y(index) = Y(index) + dt * v(2);
        end
        X(Diri) = X0(Diri);
        Y(Diri) = Y0(Diri);
        velocities(:, Diri) = 0;
        
        % plot mesh
        fig = figure(1);
        clf;
        hold on;
        quiver(X, Y, velocities(1, :)', velocities(2, :)', 'g');
        quiver(X, Y, f_elastic(1, :)', f_elastic(2, :)', 'm');
        quiver(X, Y, f_ext(1, :)', f_ext(2, :)', 'y');
        quiver(X, Y, ft(1, :)', ft(2, :)', 'k');
        triplot(T,X0,Y0,'r');
        triplot(T,X,Y,'b');
        title(t);
        xlabel('x');
        ylabel('y');
        xlim([-4 4]);
        ylim([-4, 2]);
        %axis equal;
        hold off;
        
        filename = sprintf('exp9/%d.jpg', count) ;
        saveas(fig, filename);
        
        %close(fig)
        
        t = t + dt;
        count = count + 1;
    end
    
end