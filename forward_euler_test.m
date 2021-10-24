
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