
clear all
close all
eval_u = "analytical_u_xyz";

eval_f = "analytical_f_xyz";


p.kappa=1;
p.Nx = 5;
p.Ny = 5;
p.Nz = 5;
p.hx = 1;
p.hy = 1;
p.hz = 1;

x = ones(p.Nx-1, p.Ny-1, p.Nz-1);
x = sqrt(x/numel(x));
y1 = zeros(p.Nx-1, p.Ny-1, p.Nz-1);
y2 = zeros(p.Nx-1, p.Ny-1, p.Nz-1);
y3 = zeros(p.Nx-1, p.Ny-1, p.Nz-1);

x_start = [cube2column(x);cube2column(y1);cube2column(y2);cube2column(y3)];

t_start=0;
t_stop=0.1;
visualize=1;
max_dt_FE = .005;
figure(1)
[X] = ForwardEuler(eval_f,x_start,p,eval_u,t_start,t_stop,max_dt_FE,visualize);

% PSIT = column2cube(X(1:64,t), p.Nx, p.Ny, p.Nz);
% PHITX = column2cube(X(65:128,t), p.Nx, p.Ny, p.Nz);
% PSITY = column2cube(X(129:192,t), p.Nx, p.Ny, p.Nz);
% PSITZ = column2cube(X(193:end,t), p.Nx, p.Ny, p.Nz);
n = (p.Nx-1)^3;
PSIT = X(1:n,t);
PHITX = X(n+1:2*n,t);
PSITY = X(2*n+1:3*n,t);
PSITZ = X(3*n+1:4*n,t);

[xx, yy, zz] = meshgrid(1:p.hx:p.Nx-1, 1:p.hy:p.Ny-1, 1:p.hz:p.Nz-1);
xx = cube2column(xx);
yy = cube2column(yy);
zz = cube2column(zz);

%% 
figure(2)
for t = 1:size(X, 2)
    PSIT = X(1:64,t);
    PHITX = X(65:128,t);
    PHITY = X(129:192,t);
    PHITZ = X(193:end,t);
    
    subplot(2,2,1)
    scatter3(xx,yy,zz,36,PSIT, 'filled')
    
    subplot(2,2,2)
    scatter3(xx,yy,zz,36,PHITX, 'filled')
    
    subplot(2,2,3)
    scatter3(xx,yy,zz,36,PHITY, 'filled')
    
    subplot(2,2,4)
    scatter3(xx,yy,zz,36,PHITZ, 'filled')
    
    pause(0.5)    
    sgtitle(t)
    drawnow
end
