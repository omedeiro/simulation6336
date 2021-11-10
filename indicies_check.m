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

[xx, yy, zz] = meshgrid(1:p.hx:p.hx*(p.Nx+1), 1:p.hy:p.hy*(p.Ny+1), 1:p.hz:p.hz*(p.Nz+1));


xx = cube2column(xx);
yy = cube2column(yy);
zz = cube2column(zz);


N1cube = sparse(p.m,1,0);



figure(1)
N1cube = sparse(p.m,1,0);
scatter3(xx,yy,zz,36,N1cube, 'filled', 'MarkerFaceAlpha', 0.5)
title('Nx=Ny=Nz= 10, p.m = N+1')
xlabel('x')
ylabel('y')
zlabel('z')

figure(2)
subplot(2,3,1)
N1cube = sparse(p.m,1,0);
N1cube(p.m1x) = 1;
scatter3(xx,yy,zz,36,N1cube, 'filled', 'MarkerFaceAlpha', 0.5)
title('p.m1x')
xlabel('x')
ylabel('y')
zlabel('z')

subplot(2,3,2)
N1cube = sparse(p.m,1,0);
N1cube(p.m1y) = 1;
scatter3(xx,yy,zz,36,N1cube, 'filled', 'MarkerFaceAlpha', 0.5)
title('p.m1y')
xlabel('x')
ylabel('y')
zlabel('z')

subplot(2,3,3)
N1cube = sparse(p.m,1,0);
N1cube(p.m1z) = 1;
scatter3(xx,yy,zz,36,N1cube, 'filled', 'MarkerFaceAlpha', 0.5)
title('p.m1z')
xlabel('x')
ylabel('y')
zlabel('z')

subplot(2,3,4)
N1cube = sparse(p.m,1,0);
N1cube(p.m1x_int) = 1;
scatter3(xx,yy,zz,36,N1cube, 'filled', 'MarkerFaceAlpha', 0.5)
title('p.m1x\_int')
xlabel('x')
ylabel('y')
zlabel('z')


subplot(2,3,5)
N1cube = sparse(p.m,1,0);
N1cube(p.m1y_int) = 1;
scatter3(xx,yy,zz,36,N1cube, 'filled', 'MarkerFaceAlpha', 0.5)
title('p.m1y\_int')
xlabel('x')
ylabel('y')
zlabel('z')


subplot(2,3,6)
N1cube = sparse(p.m,1,0);
N1cube(p.m1z_int) = 1;
scatter3(xx,yy,zz,36,N1cube, 'filled', 'MarkerFaceAlpha', 0.5)
title('p.m1z\_int')
xlabel('x')
ylabel('y')
zlabel('z')

%% 
figure(3)
subplot(2,3,1)
N1cube = sparse(p.m,1,0);
N1cube(p.m2x) = 1;
scatter3(xx,yy,zz,36,N1cube, 'filled', 'MarkerFaceAlpha', 0.5)
title('p.m2x')
xlabel('x')
ylabel('y')
zlabel('z')

subplot(2,3,2)
N1cube = sparse(p.m,1,0);
N1cube(p.m2y) = 1;
scatter3(xx,yy,zz,36,N1cube, 'filled', 'MarkerFaceAlpha', 0.5)
title('p.m2y')
xlabel('x')
ylabel('y')
zlabel('z')

subplot(2,3,3)
N1cube = sparse(p.m,1,0);
N1cube(p.m2z) = 1;
scatter3(xx,yy,zz,36,N1cube, 'filled', 'MarkerFaceAlpha', 0.5)
title('p.m2z')
xlabel('x')
ylabel('y')
zlabel('z')


subplot(2,3,4)
N1cube = sparse(p.m,1,0);
N1cube(p.m2x_int) = 1;
scatter3(xx,yy,zz,36,N1cube, 'filled', 'MarkerFaceAlpha', 0.5)
title('p.m2x\_int')
xlabel('x')
ylabel('y')
zlabel('z')

subplot(2,3,5)
N1cube = sparse(p.m,1,0);
N1cube(p.m2y_int) = 1;
scatter3(xx,yy,zz,36,N1cube, 'filled', 'MarkerFaceAlpha', 0.5)
title('p.m2y\_int')
xlabel('x')
ylabel('y')
zlabel('z')

subplot(2,3,6)
N1cube = sparse(p.m,1,0);
N1cube(p.m2z_int) = 1;
scatter3(xx,yy,zz,36,N1cube, 'filled', 'MarkerFaceAlpha', 0.5)
title('p.m2z\_int')
xlabel('x')
ylabel('y')
zlabel('z')

%% 

figure(4)
subplot(2,3,1)
N1cube = sparse(p.m,1,0);
N1cube(p.mNx) = 1;
scatter3(xx,yy,zz,36,N1cube, 'filled', 'MarkerFaceAlpha', 0.5)
title('p.mNx')
xlabel('x')
ylabel('y')
zlabel('z')

subplot(2,3,2)
N1cube = sparse(p.m,1,0);
N1cube(p.mNy) = 1;
scatter3(xx,yy,zz,36,N1cube, 'filled', 'MarkerFaceAlpha', 0.5)
title('p.mNy')
xlabel('x')
ylabel('y')
zlabel('z')

subplot(2,3,3)
N1cube = sparse(p.m,1,0);
N1cube(p.mNz) = 1;
scatter3(xx,yy,zz,36,N1cube, 'filled', 'MarkerFaceAlpha', 0.5)
title('p.mNz')
xlabel('x')
ylabel('y')
zlabel('z')

subplot(2,3,4)
N1cube = sparse(p.m,1,0);
N1cube(p.mNx_int) = 1;
scatter3(xx,yy,zz,36,N1cube, 'filled', 'MarkerFaceAlpha', 0.5)
title('p.mNx\_int')
xlabel('x')
ylabel('y')
zlabel('z')


subplot(2,3,5)
N1cube = sparse(p.m,1,0);
N1cube(p.mNy_int) = 1;
scatter3(xx,yy,zz,36,N1cube, 'filled', 'MarkerFaceAlpha', 0.5)
title('p.mNy\_int')
xlabel('x')
ylabel('y')
zlabel('z')


subplot(2,3,6)
N1cube = sparse(p.m,1,0);
N1cube(p.mNz_int) = 1;
scatter3(xx,yy,zz,36,N1cube, 'filled', 'MarkerFaceAlpha', 0.5)
title('p.mNz\_int')
xlabel('x')
ylabel('y')
zlabel('z')


%% 
figure(5)
subplot(2,3,1)
N1cube = sparse(p.m,1,0);
N1cube(p.mNxp1) = 1;
scatter3(xx,yy,zz,36,N1cube, 'filled', 'MarkerFaceAlpha', 0.5)
title('p.mNxp1')
xlabel('x')
ylabel('y')
zlabel('z')

subplot(2,3,2)
N1cube = sparse(p.m,1,0);
N1cube(p.mNyp1) = 1;
scatter3(xx,yy,zz,36,N1cube, 'filled', 'MarkerFaceAlpha', 0.5)
title('p.mNyp1')
xlabel('x')
ylabel('y')
zlabel('z')

subplot(2,3,3)
N1cube = sparse(p.m,1,0);
N1cube(p.mNzp1) = 1;
scatter3(xx,yy,zz,36,N1cube, 'filled', 'MarkerFaceAlpha', 0.5)
title('p.mNzp1')
xlabel('x')
ylabel('y')
zlabel('z')

subplot(2,3,4)
N1cube = sparse(p.m,1,0);
N1cube(p.mNxp1_int) = 1;
scatter3(xx,yy,zz,36,N1cube, 'filled', 'MarkerFaceAlpha', 0.5)
title('p.mNxp1\_int')
xlabel('x')
ylabel('y')
zlabel('z')


subplot(2,3,5)
N1cube = sparse(p.m,1,0);
N1cube(p.mNyp1_int) = 1;
scatter3(xx,yy,zz,36,N1cube, 'filled', 'MarkerFaceAlpha', 0.5)
title('p.mNyp1\_int')
xlabel('x')
ylabel('y')
zlabel('z')


subplot(2,3,6)
N1cube = sparse(p.m,1,0);
N1cube(p.mNzp1_int) = 1;
scatter3(xx,yy,zz,36,N1cube, 'filled', 'MarkerFaceAlpha', 0.5)
title('p.mNzp1\_int')
xlabel('x')
ylabel('y')
zlabel('z')