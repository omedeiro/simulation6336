
clear all
close all
eval_u = "analytical_u_xyz";

eval_f = "eval_f";


p.kappa = 5;
p.Nx = 100;
p.Ny = 100;
p.Nz = 100;
p.hx = 1e-2;
p.hy = 1e-2;
p.hz = 1e-2;

p = contruct_indices(p);

p.LPHIX = construct_LPHIXm(p);
p.LPHIY = construct_LPHIYm(p);
p.LPHIZ = construct_LPHIZm(p);

x = sparse(1:(p.Nx-1)*(p.Ny-1)*(p.Nz-1),1, sqrt(1/((p.Nx-1)*(p.Ny-1)*(p.Nz-1))));
y1 = sparse((p.Nx-1)*(p.Ny-1)*(p.Nz-1),1);
y2 = sparse((p.Nx-1)*(p.Ny-1)*(p.Nz-1),1);
y3 = sparse((p.Nx-1)*(p.Ny-1)*(p.Nz-1),1);

x_start = [x;y1;y2;y3];
% 
% t = 0.1;
% 
% err_J = 1e-2;
% J = eval_num_jac(x_start, eval_f, p, eval_u, t, err_J)

t_start=0;
t_stop=0.005;
max_dt_FE = .0005;

[X] = ForwardEuler(eval_f,x_start,p,eval_u,t_start,t_stop,max_dt_FE,1);


%% 
visualizeNetwork(X,p)
