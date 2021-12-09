
clear all
close all
eval_u = "eval_u";

eval_f = "eval_f";
aa = 30;
p.kappa = 5;
p.Nx = aa;
p.Ny = aa;
p.Nz = aa;
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
%{
for i = 1:(p.Nx-1)*(p.Ny-1)*(p.Nz-1)
    random = rand(1, 1);
    rphaser = rand(1, 1);
    rphasei = rand(1, 1);
    if random < 0.5
        x(i, 1) = 1 - random*(rphaser + 1i*rphasei)/(sqrt(rphaser^2 + rphasei^2));
    end
end
%}

x_start = [x;y1;y2;y3];

t_start=0;
t_stop=1;
max_dt_FE = 1E-3;
visualize = 0;

p.t_start=t_start;
p.t_stop =t_stop;
p.timestep = max_dt_FE;
p.visualizeSave = 0;

tic;
[X] = ForwardEuler(eval_f,x_start,p,eval_u,t_start,t_stop,max_dt_FE,visualize);
toc

p.frames = 5;
p.visualizeSave=1;
visualizeNetworkX(X, p)

%{
%% 
OM=0;
if OM == 1
[Bx,By,Bz] = eval_Bfield(X,p);

[yy, xx, zz] = meshgrid(1:p.hx:p.hx*(p.Nx-1), 1:p.hy:p.hy*(p.Ny-1), 1:p.hz:p.hz*(p.Nz-1));

[yy2, xx2, zz2] = meshgrid(p.hx:p.hx:p.hx*(p.Nx-2), p.hy:p.hy:p.hy*(p.Ny-2), p.hz:p.hz:p.hz*(p.Nz-2));
xx = cube2column(xx);
yy = cube2column(yy);
zz = cube2column(zz);
xx2 = cube2column(xx2);
yy2 = cube2column(yy2);
zz2 = cube2column(zz2);


n = (p.Nx-1)*(p.Ny-1)*(p.Nz-1);
n2 = (p.Nx-2)*(p.Ny-2)*(p.Nz-2);
figure
for i = 1:p.Nx*p.Ny+p.Nz
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

end
%}
