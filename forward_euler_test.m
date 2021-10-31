
clear all
close all
eval_u = "analytical_u_xyz";

eval_f = "eval_f";


p.kappa=1;
p.Nx = 20;
p.Ny = 20;
p.Nz = 20;
p.hx = 1;
p.hy = 1;
p.hz = 1;

x = sparse((p.Nx-1)*(p.Ny-1)*(p.Nz-1),1, sqrt(1/((p.Nx-1)*(p.Ny-1)*(p.Nz-1))));
y1 = sparse((p.Nx-1)*(p.Ny-1)*(p.Nz-1),1);
y2 = sparse((p.Nx-1)*(p.Ny-1)*(p.Nz-1),1);
y3 = sparse((p.Nx-1)*(p.Ny-1)*(p.Nz-1),1);


x_start = [x;y1;y2;y3];

t_start=0;
t_stop=0.010;
visualize=0;
max_dt_FE = .001;
% figure(1)
[X] = ForwardEuler(eval_f,x_start,p,eval_u,t_start,t_stop,max_dt_FE,visualize);


%% 


[xx, yy, zz] = meshgrid(1:p.hx:p.Nx-1, 1:p.hy:p.Ny-1, 1:p.hz:p.Nz-1);
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
