clearvars;
close all;
clc;

% load('fvm_data.mat');  % Load mesh
% fd= @(p) drectangle(p,-3,3,-3,3);
% es = 0.2;
% [p, T] = distmesh2d(fd,@huniform,es,[-3,-3;3,3],[-3,-3; -3,3; 3,-3; 3,3]);
fd= @(p) drectangle(p,-3,3,-1,1);
es = 0.5;
[p, T] = distmesh2d(fd,@huniform,es,[-3,-3;3,3],[-3,-1; -3,1; 3,-1; 3,1]);
X      = p ( : , 1 ) ;
Y      = p ( : , 2 ) ;
cntV   = size(X,1);
V      = (1:cntV)';
cntT   = size(T,1);

dx1 = X( T(:,3) ) - X( T(:,2) );
dy1 = Y( T(:,3) ) - Y( T(:,2) );
dx2 = X( T(:,1) ) - X( T(:,2) );
dy2 = Y( T(:,1) ) - Y( T(:,2) );
A =  (dx1.*dy2 - dx2.*dy1 ) ./ 2; % Compute triangle areas

%%
% [TR, CVs, Bmask] = create_control_volumes(T, X, Y);

% figure(1);
% clf;
% draw_control_volumes(TR, CVs, 0.05);
% axis equal;
% title('Computational mesh')

% [A, b] = matrix_assembly(T, X, Y, CVs);

%%
E  = 69e9;    % Young modulus
nu = 0.3;    % Poisson ration

lambda = 1000;
mu = 1000;

% lambda = (E * nu) / ((1 + nu) * (1 - 2 * nu));
% mu = E / (2 * (1 + nu)); 

% mu = 0.0006; % rubber
% mu = 79.3; % steel
% nu = 0.4999; % rubber
% nu = 0.3; % steel
% lambda = (2 * mu * nu) / (1 - 2 * nu);

% E = 0.01e9; % rubber
% nu = 0.48; % rubber 
% rho = 1050; % rubber
% 
% mu = E / (2 * (1 + nu));
% lambda = (E * nu) / ((1 + nu) * (1 - 2 * nu));

x = X;
y = Y;
dt = 0.001;
% dt = 1 / E;
% m = rho;
m= 1;
v = zeros(cntV,2);
% fe = zeros(cntV,2);
% tforce(100:2000) = 0;
% iters = 8;
iters = 2000;
% iters = 1 / dt;
t_force = -ones(1,iters) * 20;%1000000;
fixedb = find(X > -2.99);
lfb = length(fixedb);


TR2 = triangulation(T,X,Y);
Ts = vertexAttachments(TR2);
es = edges(TR2);
es1 = es(:,1);
es2 = es(:,2);
b = es(X(es1) > 2.99 & X(es2) > 2.99,:);

d = @(x,y) sqrt(sum((x - y).^2));

D0 = zeros(2,2,cntT);
ni = zeros(2,cntV,6);
a_e = zeros(cntV,1);

for t = 1:cntT
    i = T(t,1);
    j = T(t,2);
    k = T(t,3);
    D0(:,:,t) = [X(j) - X(i) X(k) - X(i) ; Y(j) - Y(i) Y(k) - Y(i)]^-1;
end

for i = 1:cntV
    ts = Ts{i};
    for e = 1:length(ts)
        t = T(ts(e),:);
        [~,j,k] = find_vertex_order(i,t(1),t(2),t(3));
        ni(:,i,e) = (-1 * -[Y(i) - Y(j) ; X(j) - X(i)] -1 * [Y(i) - Y(k) ; X(k) - X(i)]) / 2;
        a_e(i) = a_e(i) + A(ts(e));
    end
end

wb = waitbar(0,'');

for step = 1:iters
    fe = zeros(2,cntV);
    Pe_a = zeros(2,2,cntT);
    
    for t = 1:cntT
        i = T(t,1);
        j = T(t,2);
        k = T(t,3);
        D = [x(j) - x(i) x(k) - x(i) ; y(j) - y(i) y(k) - y(i)];
%         Fe = D * [X(j) - X(i) X(k) - X(i) ; Y(j) - Y(i) Y(k) - Y(i)]^-1;
        Fe = D * D0(:,:,t);
        Ee = (Fe' * Fe - [1 0 ; 0 1]) / 2;
        Se = lambda * trace(Ee) * [1 0 ; 0 1] + 2 * mu * Ee;
        Pe_a(:,:,t) = Fe * Se;
    end

    % CVs = cell(1,cntV);
    for i = 1:cntV
        ts = Ts{i};
%         a = zeros(3,length(ts));
        %     fe(i,:) = [0; 0];
%         a_e(i) = 0;
        fei = [0 ; 0];
%         a_ei = 0;
%         nie = ni(:,:,i);
        for e = 1:length(ts)
%             t = T(ts(e),:);
%             [~,j,k] = find_vertex_order(i,t(1),t(2),t(3));
%             D = [x(j) - x(i) x(k) - x(i) ; y(j) - y(i) y(k) - y(i)];
%             Fe = D * [X(j) - X(i) X(k) - X(i) ; Y(j) - Y(i) Y(k) - Y(i)]^-1;
%             Ee = (Fe' * Fe - eye(2)) / 2;
%             Se = lambda * trace(Ee) * eye(2) + 2 * mu * Ee;
%             Pe = Fe * Se;
            Pe = Pe_a(:,:,ts(e));
%             ni(:,i,e) = (-1 * -[Y(i) - Y(j) ; X(j) - X(i)] -1 * [Y(i) - Y(k) ; X(k) - X(i)]) / 2;
            fei = fei + Pe * ni(:,i,e);%ni(:,e,i);
%             a_ei = a_ei + A(ts(e));
            %         a(:,e) = [i, j, k];
        end
        fe(:,i) = fei;
%         a_e(i) = a_ei;
        %     CVs{i} = a
    end
    
    for e = 1:(size(b,1))
        xi = X(b(e,1));
        xj = X(b(e,2));
        yi = Y(b(e,1));
        yj = Y(b(e,2));
        fl = [0 t_force(step)]';
        l_half = d([xi yi],[xj yj])/2;
        fe(:,b(e,1)) = fe(:,b(e,1)) + fl * l_half;
        fe(:,b(e,2)) = fe(:,b(e,2)) + fl * l_half;
        %     f_t(1:2,e) = fl * l_half;
    end
    
    f_total = fe';
    
%     for i = 1:cntV
%         v(i,:) = v(i,:) + (dt / (m * a_e(i) / 3)) * f_total(i,:);
%         v(X < -3 + 0.01,:) = 0;
%         x(i) = x(i) + dt * v(i,1);
%         y(i) = y(i) + dt * v(i,2);
%     end
    v = v + (dt ./ (m * a_e / 3)) .* f_total;
    v(X < -3 + 0.01,:) = 0;
    x = x + dt .* v(:,1);
    y = y + dt .* v(:,2);
    x_t(step,:) = x;
    y_t(step,:) = y;
%     z = lsqnonlin(@(z) f(x,y,v,z,dt,m,cntV,cntT,T,Ts,D0,ni,b,X,Y,fixedb,lfb),zeros(1,lfb*4));
%     x(fixedb) = z(1:lfb);
%     y(fixedb) = z(lfb+1:2*lfb);
%     v(fixedb,1) = z(2*lfb+1:3*lfb);
%     v(fixedb,2) = z(3*lfb+1:end);
    x_t(step,:) = x;
    y_t(step,:) = y;
%     triplot(T,x,y);
%     axis equal;
%     drawnow;
    waitbar(step/iters,wb,'');
end


close(wb);

for step = 1:16:iters
    triplot(T,x_t(step,:),y_t(step,:));
%     axis equal;
    axis([-4 4 -4 4]);
    drawnow;
end

% triplot(T,X,Y,'b');
% hold on;
% triplot(T,x,y,'r');
% hold off;
% title('The computational mesh');
% xlabel('x');
% ylabel('y');
% axis equal;

function r = f(x,y,v,z,dt,m,cntV,cntT,T,Ts,D0,ni,b,X,Y,fixedb,lfb)
lambda = 1000;
mu = 1000;
tf = -20;
N = lfb;
xt = [x,y];
xt(fixedb,1) = z(1:N);
xt(fixedb,2) = z(N+1:N+N);
vt = zeros(cntV,2);
vt(fixedb,1) = z(2*N+1:3*N);
vt(fixedb,2) = z(3*N+1:end);
% vt(:,1) = z(2*N+1:3*N);
% vt(:,2) = z(3*N+1:4*N);
% fx = zeros(2,N)';
    fe = zeros(2,cntV);
    Pe_a = zeros(2,2,cntT);
    
    for t = 1:cntT
        i = T(t,1);
        j = T(t,2);
        k = T(t,3);
        D = [x(j) - x(i) x(k) - x(i) ; y(j) - y(i) y(k) - y(i)];
%         Fe = D * [X(j) - X(i) X(k) - X(i) ; Y(j) - Y(i) Y(k) - Y(i)]^-1;
        Fe = D * D0(:,:,t);
        Ee = (Fe' * Fe - [1 0 ; 0 1]) / 2;
        Se = lambda * trace(Ee) * [1 0 ; 0 1] + 2 * mu * Ee;
        Pe_a(:,:,t) = Fe * Se;
    end

    % CVs = cell(1,cntV);
    for i = 1:cntV
        ts = Ts{i};
%         a = zeros(3,length(ts));
        %     fe(i,:) = [0; 0];
%         a_e(i) = 0;
        fei = [0 ; 0];
%         a_ei = 0;
%         nie = ni(:,:,i);
        for e = 1:length(ts)
%             t = T(ts(e),:);
%             [~,j,k] = find_vertex_order(i,t(1),t(2),t(3));
%             D = [x(j) - x(i) x(k) - x(i) ; y(j) - y(i) y(k) - y(i)];
%             Fe = D * [X(j) - X(i) X(k) - X(i) ; Y(j) - Y(i) Y(k) - Y(i)]^-1;
%             Ee = (Fe' * Fe - eye(2)) / 2;
%             Se = lambda * trace(Ee) * eye(2) + 2 * mu * Ee;
%             Pe = Fe * Se;
            Pe = Pe_a(:,:,ts(e));
%             ni(:,i,e) = (-1 * -[Y(i) - Y(j) ; X(j) - X(i)] -1 * [Y(i) - Y(k) ; X(k) - X(i)]) / 2;
            fei = fei + Pe * ni(:,i,e);%ni(:,e,i);
%             a_ei = a_ei + A(ts(e));
            %         a(:,e) = [i, j, k];
        end
        fe(:,i) = fei;
%         a_e(i) = a_ei;
        %     CVs{i} = a
    end
    
    for e = 1:(size(b,1))
        xi = X(b(e,1));
        xj = X(b(e,2));
        yi = Y(b(e,1));
        yj = Y(b(e,2));
        fl = [0 tf]';
        l_half = sum(([xi yi] - [xj yj]).^2)/2;
        fe(:,b(e,1)) = fe(:,b(e,1)) + fl * l_half;
        fe(:,b(e,2)) = fe(:,b(e,2)) + fl * l_half;
        %     f_t(1:2,e) = fl * l_half;
    end
z = [vt - v - (dt/m) * fe', xt - [x y] - dt * vt];
z = z(:);
r = z;
% r = sum(z.^2);
end
