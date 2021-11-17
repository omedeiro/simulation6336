
clear all
close all

eval_f = "eval_f";
eval_u = "eval_u";

p.kappa = 5;
p.Nx = 10;
p.Ny = 10;
p.Nz = 10;
p.hx = 1;
p.hy = 1;
p.hz = 1;

p.magBx = 0;
p.magBy = 0;
p.magBz = 0;

p.periodic_x = 0;
p.periodic_y = 0;
p.periodic_z = 0;

p = contruct_indices(p);

p.LPHIX = construct_LPHIXm(p);
p.LPHIY = construct_LPHIYm(p);
p.LPHIZ = construct_LPHIZm(p);

x = sparse(1:(p.Nx-1)*(p.Ny-1)*(p.Nz-1),1, 1);
y1 = sparse((p.Nx-1)*(p.Ny-1)*(p.Nz-1),1);
y2 = sparse((p.Nx-1)*(p.Ny-1)*(p.Nz-1),1);
y3 = sparse((p.Nx-1)*(p.Ny-1)*(p.Nz-1),1);

x_start = [x;y1;y2;y3];

p.t_start=0;
p.t_stop=.1;
p.timestep = 1e-2;
visualize = 1;
p.visualizeSave = 0;
[X] = Trapezoidal(eval_f,x_start,p,eval_u,p.t_start,p.t_stop,p.timestep,visualize);



