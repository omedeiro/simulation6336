%%%%%%% SENSITIVITY %%%%%%%

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

X0 = [cube2column(x);cube2column(y1);cube2column(y2);cube2column(y3)];

%define F
F = @(X) analytical_f_xyz(X, Bx, hx, hy, hz, kappa, Nx, Ny, Nz);

tstart = tic;
J = eval_num_jac(X0, F);
telapsed = toc(tstart);

% imagesc(abs(J))
X = J\(-F(X0));

N = (Nx-1)*(Ny-1)*(Nz-1);
c_x = [ones(N,1)/N; zeros(3*N,1)];
c_y1 = [zeros(N,1); ones(N,1)/N; zeros(2*N,1)];
c_y2 = [zeros(2*N,1); ones(N,1)/N; zeros(N,1)];
c_y3 = [zeros(3*N,1); ones(N,1)/N];
c_y = [zeros(N,1); ones(3*N,1)/(3*N)];            % average of all y
    
cube_FrontFace = zeros(Nx-1, Ny-1, Nz-1);
cube_FrontFace(1,:,:) = ones(Ny-1, Nz-1)/(2*(Nz-1)*(Ny-1));
c_FrontFace = [zeros(2*N,1); cube2column(cube_FrontFace); cube2column(cube_FrontFace)]; % average of y2 and y3 in the front face along x

C = [c_x c_y1 c_y2 c_y3 c_y c_FrontFace];

Y = C'*X;
p = kappa;
dp_p = [-99 -80 -70 -60 -50 -40 -30 -20 -10 -5 1 5 10 20 30 40 50 60 70 80 90 100]; % dp/p in %
dY_dp = zeros(size(C,2), numel(dp_p));
dY_Y = zeros(size(C,2), numel(dp_p));
dJ_J = zeros(1, numel(dp_p));

for i = 1 : numel(dp_p) 
    eps = (dp_p(i)/100)*p;
    F_eps = @(X) analytical_f_xyz(X, Bx, hx, hy, hz, p+eps, Nx, Ny, Nz);
    J_eps = eval_num_jac(X0, F_eps);
    X_eps = J_eps\(-F_eps(X0));

    dY_dp(:,i) = C'*(X_eps - X)/eps;

    % percentage variation
    dY_Y(:,i) = 100*eps*abs(dY_dp(:,i))./abs(Y);
    dJ_J(i) = 100*norm(abs(J_eps - J))/norm(abs(J));
end

plot(dp_p, abs(dY_Y'), 'LineWidth', 1)
set(gca, 'YScale', 'log')
legend('|y_x|', '|y_{y1}|', '|y_{y2}|', '|y_{y3}|', '|y_{y}|', '|y_{FrontFace}|')
ylabel('% variation of |y|')
xlabel('% variation of k')
grid on

Y
dY_dp 
dY_Y 
dJ_J

%%%%%%%%%%%adjoint not convenient because 6 Ys and only 1 param