Nx = 3;
Ny = 3;
Nz = 3;
x = ones(Nx-1, Ny-1, Nz-1);
y1 = zeros(Nx-1, Ny-1, Nz-1);
y2 = zeros(Nx-1, Ny-1, Nz-1);
y3 = zeros(Nx-1, Ny-1, Nz-1);
hx = 1;
hy = 1;
hz = 1;
kappa=1;
Bx=1;

U0=Bx;
X0 = [cube2column(x);cube2column(y1);cube2column(y2);cube2column(y3)];

%define F
F = @(X) analytical_f_xyz(X, Bx, hx, hy, hz, kappa, Nx, Ny, Nz);
FU = @(U) analytical_f_xyz(X0, U, hx, hy, hz, kappa, Nx, Ny, Nz);

% JnumericU = eval_num_jacU(U0, FU);

tstart = tic;
Jnumeric = eval_num_jac(X0, F, 1e-6);
JnumericU = eval_num_jac(U0, FU, 1e-6);


imagesc(abs(Jnumeric))



% Janalytic = eval_analytical_jac_xyz(x, y1, y2, y3, hx, hy, hz, kappa, Nx-1);
% 