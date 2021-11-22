%Linearization


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
