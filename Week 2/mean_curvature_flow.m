function phi = mean_curvature_flow(phi,T)
% MEAN_CURVATURE_FLOW - Mean Curvature Flow
%   psi = curvatureflow(phi,T)
% INPUT:
%      phi - Initial signed distance map
%        T - Total running time
% RESULT:
%      psi - Resulting level set surface
%
    dx = 1;
    dy = dx;
    h = min(dx, dy);
    kappa = 1 / max(dx, dy);
    dt = h/2*kappa;
    l = size(phi, 1);
    
    % compute Dx, Dy, Dxx, Dyy, Dxy
    k = zeros(l);
    for t = 0:dt:T
         % function to fill phi with ghost nodes
        phi = ghost_nodes(phi);
        for i = 2:l + 1
            for j = 2:l + 1
                Dx = (phi(i + 1, j) - phi(i - 1, j)) / (2*dx);
                Dy = (phi(i, j + 1) - phi(i, j - 1)) / (2*dy);
                Dxx = (phi(i + 1, j) - 2*phi(i, j) + phi(i - 1, j)) / (dx^2);
                Dyy = (phi(i, j + 1) - 2*phi(i, j) + phi(i, j - 1)) / (dy^2);
                Dxy = (phi(i + 1, j + 1) - phi(i + 1, j - 1) - phi(i - 1, j + 1) + phi(i - 1, j - 1)) / (4*dx*dy);

                temp = (Dx)^2 * Dyy + (Dy)^2 * Dxx - 2*Dxy * Dx * Dy;
                g = sqrt((Dx)^2 + (Dy)^2);

                % numerical remedies
                if g <= 0.5
                    g = 1;
                end
                temp = temp/g^3;

                k(i - 1, j - 1) = max(-kappa, min(temp, kappa));
            end
        end
        
        % compute new phi 
        phi = phi(2:l + 1, 2:l + 1) + dt * k; 
        
        imagesc(phi); colormap(gray); axis image; colorbar
        hold on; title(t);contour(phi,[0,0],'r'); hold off
        drawnow;
    end

end
