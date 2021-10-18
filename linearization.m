%Linearization

Nx = 3;
Ny = 3;
Nz = 3;
x = ones(Nx-1, Ny-1, Nz-1);
x = sqrt(x/numel(x));
y1 = ones(Nx-1, Ny-1, Nz-1);
y2 = ones(Nx-1, Ny-1, Nz-1);
y3 = ones(Nx-1, Ny-1, Nz-1);
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



J0F = eval_num_jac(X0,F,1e-5);
A = J0F;

J0U = eval_num_jac(U0,FU, 1e-5);
K0 = F(X0) - J0F*X0 - J0U*U0;
B = [K0 J0U];

xeval = J0F\-F(X0);
linf = A*X0 + B*[1;U0];
f = F(xeval)
imagesc(abs(f-linf))
colorbar
