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
Bx=0;


X0 = [cube2column(x);cube2column(y1);cube2column(y2);cube2column(y3)];

%define F
F = @(X) analytical_f_xyz(X, Bx, hx, hy, hz, kappa, Nx, Ny, Nz);

tstart = tic;
Jnumeric = eval_num_jac(X0, F);
telapsed = toc(tstart);


Janalytic = eval_analytical_jac_xyz(x, y1, y2, y3, hx, hy, hz, kappa, Nx-1);

function J = eval_num_jac(X0, F)
    eps_Im = .01;
    eps_Re = .01;
    S_Im = 2;
    S_Re = 2;
    J = zeros(size(F(X0),1), numel(X0));
    Jp = ones(size(F(X0),1), numel(X0));
    err = 1e-6;
    while any(abs((J - Jp)) > err, 'all') 
        Jp = J;
        for k = 1 : size(J,2) % loop columns
            dx = zeros(size(X0,1), 1);
            dx(k) = eps_Re + 1i*eps_Im;
            J(:, k) = (F(X0 + dx) - F(X0))/dx(k);
        end
        eps_Re = eps_Re/S_Re;
        eps_Im = eps_Im/S_Im;
        disp(max(max(abs((J - Jp)))))
    end
end
