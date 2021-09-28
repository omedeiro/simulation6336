N = 3;
psi = [1; 0; 0];
phix = [0; 0; 0];

hx = 1;

state_0 = cat(1, psi, phix);

%define F
F = @(x) analytical_f(x, hx, N);

Jnumeric = eval_num_jac(state_0, F)
Janalytic = eval_analytical_jac(state_0, hx, N)

imagesc(abs(Janalytic-Jnumeric))
