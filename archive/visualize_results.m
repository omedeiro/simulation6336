Nx = 4;
Ny = 4;
Nz = 4;

x = ones(Nx-1, Ny-1, Nz-1);
y1 = ones(Nx-1, Ny-1, Nz-1);
y2 = ones(Nx-1, Ny-1, Nz-1);
y3 = ones(Nx-1, Ny-1, Nz-1);
hx = 1;
hy = 1;
hz = 1;
kappa=1;
Bx=1;

[xx,yy,zz] = meshgrid(1:Nx, 1:Ny, 1:Nz)
px = [xx(:);xx(:);xx(:)];
py = [yy(:);yy(:);yy(:)];
pz = [zz(:);zz(:);zz(:)];

X0 = [cube2column(x);cube2column(y1);cube2column(y2);cube2column(y3)];

%define F
F = @(X) analytical_f_xyz(X, Bx, hx, hy, hz, kappa, Nx, Ny, Nz);

tstart = tic;
J = eval_num_jac(X0, F, 1e-5);
telapsed = toc(tstart);

% imagesc(abs(J))
X = J\(-F(X0));