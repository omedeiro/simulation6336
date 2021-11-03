
clear all
close all
eval_u = "analytical_u_xyz";

eval_f = "eval_f";


p.kappa=1;
p.Nx = 30;
p.Ny = 30;
p.Nz = 30;
p.hx = 1e-9;
p.hy = 1e-9;
p.hz = 1e-9;


h = 1;
for k = 2 : p.Nz
    for j = 2 : p.Ny
        for i = 2 : p.Nx
            p.M2(h) = i + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
            p.m2(h) = i-1 + (p.Nx-1)*(j-2)+(p.Nx-1)*(p.Ny-1)*(k-2);
            h = h + 1;
        end
    end
end
% p.M2 = sub2ind([(p.Nx+1) (p.Nx+1) (p.Ny+1)], [2:p.Nx], [2:p.Ny], [2:p.Nz]);
% p.m2 = sub2ind([(p.Nx-1) (p.Nx-1) (p.Ny-1)], [1:p.Nx-1], [1:p.Ny-1], [1:p.Nz-1]);

h = 1;
for k = 1 : p.Nz+1
    for j = 1 : p.Ny+1
        for i = 1 : p.Nx+1
            p.m(h) = i + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
            h = h + 1;
        end
    end
end
% 
% p.M2 = sub2ind([(p.Nx+1) (p.Nx+1) (p.Ny+1)], [1:p.Nx+1], [1:p.Ny+1], [1:p.Nz+1]);

% x
h = 1;
h_int = 1;
for k = 1 : p.Nz+1
    for j = 1 : p.Ny+1
        p.m1x(h) = 1 + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
        p.m2x(h) = 2 + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
        p.mNx(h) = p.Nx + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
        p.mNxp1(h) = p.Nx+1 + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
        h = h + 1;
        if k ~= 1 && k ~= p.Nz+1 && j ~= 1 && j ~= p.Ny+1
            p.m1x_int(h_int) = 1 + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
            p.m2x_int(h_int) = 2 + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
            p.mNx_int(h_int) = p.Nx + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
            p.mNxp1_int(h_int) = p.Nx+1 + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
            h_int = h_int + 1;
        end           
    end
end


% y
h = 1;
h_int = 1;
for k = 1 : p.Nz+1
    for i = 1 : p.Nx+1 
        p.m1y(h) = i + (p.Nx+1)*(p.Ny+1)*(k-1);
        p.m2y(h) = i + (p.Nx+1)+(p.Nx+1)*(p.Ny+1)*(k-1);
        p.mNy(h) = i + (p.Nx+1)*(p.Ny-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
        p.mNyp1(h) = i + (p.Nx+1)*p.Ny+(p.Nx+1)*(p.Ny+1)*(k-1);
        h = h + 1;
        if k ~= 1 && k ~= p.Nz+1 && i ~= 1 && i ~= p.Nx+1
            p.m1y_int(h_int) = i + (p.Nx+1)*(p.Ny+1)*(k-1);
            p.m2y_int(h_int) = i + (p.Nx+1)+(p.Nx+1)*(p.Ny+1)*(k-1);
            p.mNy_int(h_int) = i + (p.Nx+1)*(p.Ny-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
            p.mNyp1_int(h_int) = i + (p.Nx+1)*p.Ny+(p.Nx+1)*(p.Ny+1)*(k-1);
            h_int = h_int + 1;
        end           
    end
end

% z
h = 1;
h_int = 1;
for j = 1 : p.Ny+1
    for i = 1 : p.Nx+1
        p.m1z(h) = i + (p.Nx+1)*(j-1);
        p.m2z(h) = i + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1);
        p.mNz(h) = i + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(p.Nz-1);
        p.mNzp1(h) = i + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*p.Nz;
        h = h + 1;
        if j ~= 1 && j ~= p.Ny+1 && i ~= 1 && i ~= p.Nx+1
            p.m1z_int(h_int) = i + (p.Nx+1)*(j-1);
            p.m2z_int(h_int) = i + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1);
            p.mNz_int(h_int) = i + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(p.Nz-1);
            p.mNzp1_int(h_int) = i + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*p.Nz;
            h_int = h_int + 1;
        end
    end
end

x = sparse((p.Nx-1)*(p.Ny-1)*(p.Nz-1),1, sqrt(1/((p.Nx-1)*(p.Ny-1)*(p.Nz-1))));
y1 = sparse((p.Nx-1)*(p.Ny-1)*(p.Nz-1),1);
y2 = sparse((p.Nx-1)*(p.Ny-1)*(p.Nz-1),1);
y3 = sparse((p.Nx-1)*(p.Ny-1)*(p.Nz-1),1);

x_start = [x;y1;y2;y3];

t_start=0;
t_stop=0.003;
visualize=0;
max_dt_FE = .001;
% figure(1)
[X] = ForwardEuler(eval_f,x_start,p,eval_u,t_start,t_stop,max_dt_FE,visualize);



%% 


[xx, yy, zz] = meshgrid(0:p.hx:p.hx*(p.Nx-2), 0:p.hy:p.hy*(p.Ny-2), 0:p.hz:p.hz*(p.Nz-2));
xx = cube2column(xx);
yy = cube2column(yy);
zz = cube2column(zz);
n = (p.Nx-1)*(p.Ny-1)*(p.Nz-1);

for t = 1:size(X, 2)
    figure(2)
    PSIT = X(1:n,t);
    PHITX = X(n+1:2*n,t);
    PHITY = X(2*n+1:3*n,t);
    PHITZ = X(3*n+1:4*n,t);

    subplot(2,2,1)
    scatter3(xx,yy,zz,36,real(PSIT), 'filled', 'MarkerFaceAlpha', 0.5)
    colorbar
    axis equal
    xlabel('x')
    ylabel('y')
    zlabel('z')
    
    subplot(2,2,2)
    scatter3(xx,yy,zz,36,real(PHITX), 'filled', 'MarkerFaceAlpha', 0.5)
    colorbar
    axis equal
    xlabel('x')
    ylabel('y')
    zlabel('z')
    
    subplot(2,2,3)
    scatter3(xx,yy,zz,36,real(PHITY), 'filled', 'MarkerFaceAlpha', 0.5)
    colorbar
    axis equal
    xlabel('x')
    ylabel('y')
    zlabel('z')
    
    subplot(2,2,4)
    scatter3(xx,yy,zz,36,abs(PHITZ), 'filled', 'MarkerFaceAlpha', 0.5)
    colorbar
    axis equal
    xlabel('x')
    ylabel('y')
    zlabel('z')
    
    pause(0.5)    
    sgtitle(t)
    drawnow
    
    figure(3)
    PSIT = X(1:n,t);
    PHITX = X(n+1:2*n,t);
    PHITY = X(2*n+1:3*n,t);
    PHITZ = X(3*n+1:4*n,t);

    subplot(2,2,1)
    scatter3(xx,yy,zz,36,imag(PSIT), 'filled', 'MarkerFaceAlpha', 0.5)
    colorbar
    axis equal
    xlabel('x')
    ylabel('y')
    zlabel('z')
    
    subplot(2,2,2)
    scatter3(xx,yy,zz,36,imag(PHITX), 'filled', 'MarkerFaceAlpha', 0.5)
    colorbar
    axis equal
    xlabel('x')
    ylabel('y')
    zlabel('z')
    
    subplot(2,2,3)
    scatter3(xx,yy,zz,36,real(PHITY), 'filled', 'MarkerFaceAlpha', 0.5)
    colorbar
    axis equal
    xlabel('x')
    ylabel('y')
    zlabel('z')
    
    subplot(2,2,4)
    scatter3(xx,yy,zz,36,real(PHITZ), 'filled', 'MarkerFaceAlpha', 0.5)
    colorbar
    axis equal
    xlabel('x')
    ylabel('y')
    zlabel('z')
    
    pause(0.05)    
    sgtitle(t)
    drawnow
end
