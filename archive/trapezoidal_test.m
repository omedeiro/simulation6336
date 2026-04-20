
clear all
close all

% eval_f = "eval_f_B";
% eval_u = "analytical_u_xyz";

eval_f = "eval_f";
eval_u = "eval_u";


p.kappa = 5;
p.Nx = 60;
p.Ny = 60;
p.Nz = 3;
p.hx = 1;
p.hy = 1;
p.hz = 1;

p.magBx = 0;
p.magBy = 0;
p.magBz = 0.6;
p.appliedBz = 0;
p.periodic_x = 0;
p.periodic_y = 0;
p.periodic_z = 0;
p.slice=1;
p = contruct_indices(p);

p.LPHIX = construct_LPHIXm(p);
p.LPHIY = construct_LPHIYm(p);
p.LPHIZ = construct_LPHIZm(p);

if p.Nz > 1
    x = sparse(1:(p.Nx-1)*(p.Ny-1)*(p.Nz-1),1, 1);
    y1 = sparse((p.Nx-1)*(p.Ny-1)*(p.Nz-1),1);
    y2 = sparse((p.Nx-1)*(p.Ny-1)*(p.Nz-1),1);
    y3 = sparse((p.Nx-1)*(p.Ny-1)*(p.Nz-1),1);
else
    x = sparse(1:(p.Nx-1)*(p.Ny-1),1, 1);
    y1 = sparse((p.Nx-1)*(p.Ny-1),1);
    y2 = sparse((p.Nx-1)*(p.Ny-1),1);
    y3 = [];
end

% initialize as quenched normal state

% for i = 1:(p.Nx-1)*(p.Ny-1)*(p.Nz-1)
%     random = rand(1, 1);
%     rphaser = rand(1, 1);
%     rphasei = rand(1, 1);
%     if random < 0.5
%         x(i, 1) = 1 - random*(rphaser + 1i*rphasei)/(sqrt(rphaser^2 + rphasei^2));
%     end
% end

% initialize w/vortex
%{
offset = 7;
for i=1:8
    for j = 1:8
        x(index_map(offset+i, offset+j, 1, p)) = 0;
    end
end
%}


%%
x_start = [x;y1;y2;y3];
p.linearize = 0;
p.t_start=0;
p.t_stop=100;
p.timestep = 1e-1;
visualize = 0;
p.visualizeSave = 0;
[X,p] = Trapezoidal(eval_f,x_start,p,eval_u,p.t_start,p.t_stop,p.timestep,visualize);


%% 
p.frames = 250;
p.visualizeSave=1;
visualizeNetworkX2d(X,p)
% save('Xom.mat','X')