
clear all
close all
eval_u = "analytical_u_xyz";

eval_f = "eval_f";

p.kappa = 5;
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

% x = sparse(1:(p.Nx-1)*(p.Ny-1)*(p.Nz-1),1, sqrt(1/((p.Nx-1)*(p.Ny-1)*(p.Nz-1))));
x = sparse(1:(p.Nx-1)*(p.Ny-1)*(p.Nz-1),1, 1);
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
t_stop=.005;
max_dt_FE = .0005;
visualize = 0;

[X] = ForwardEuler(eval_f,x_start,p,eval_u,t_start,t_stop,max_dt_FE,visualize);

%%
visualizeNetworkFace(X,p)

%% 
[Bx,By,Bz] = eval_Bfield(X,p);

[xx, yy, zz] = meshgrid(1:p.hx:p.hx*(p.Nx-1), 1:p.hy:p.hy*(p.Ny-1), 1:p.hz:p.hz*(p.Nz-1));

[xx2, yy2, zz2] = meshgrid(p.hx:p.hx:p.hx*(p.Nx-2), p.hy:p.hy:p.hy*(p.Ny-2), p.hz:p.hz:p.hz*(p.Nz-2));
xx = cube2column(xx);
yy = cube2column(yy);
zz = cube2column(zz);
xx2 = cube2column(xx2);
yy2 = cube2column(yy2);
zz2 = cube2column(zz2);


n = (p.Nx-1)*(p.Ny-1)*(p.Nz-1);
n2 = (p.Nx-2)*(p.Ny-2)*(p.Nz-2);
figure
for i = 1:length(p.m)
    scatter3(xx2(i),yy2(i),zz2(i),36,abs(Bx(i,1)), 'filled', 'MarkerFaceAlpha', 0.5)
    colorbar
%     caxis([0 max(max(Bx))]);
    axis equal
    title('Bx')
    xlabel('x')
    ylabel('y')
    zlabel('z')
    drawnow
    hold on
    pause(0.05)
end
