clear all
close all
eval_u = "analytical_u_xyz";

eval_f = "eval_f";


p.kappa = 1;
p.Nx = 10;
p.Ny = 10;
p.Nz = 10;
p.hx = 1;
p.hy = 1;
p.hz = 1;

p = contruct_indices(p);

p.LPHIX = construct_LPHIXm(p);
p.LPHIY = construct_LPHIYm(p);
p.LPHIZ = construct_LPHIZm(p);

x = sparse(1:(p.Nx-1)*(p.Ny-1)*(p.Nz-1),1, 1);
y1 = sparse((p.Nx-1)*(p.Ny-1)*(p.Nz-1),1);
y2 = sparse((p.Nx-1)*(p.Ny-1)*(p.Nz-1),1);
y3 = sparse((p.Nx-1)*(p.Ny-1)*(p.Nz-1),1);

x_start = [x;y1;y2;y3];



u = feval(eval_u, 0);
p.errf = 1E-3;
p.errDeltax = 1E-3;
p.relDeltax = 1E-3;
p.MaxIter = 20;
p.tolrGCR = 1E-3;
p.epsMF = 1E-3;


t_start=0;
t_stop=.005;
max_dt_FE = .0005;

[X] = Trapezoidal(eval_f,x_start,p,eval_u,t_start,t_stop,max_dt_FE,1);


%% 
% visualizeNetwork(1, X,p)
