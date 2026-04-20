%Linearization

clear all
close all

% eval_f = "eval_f_B";
% eval_u = "analytical_u_xyz";

p.eval_f = "eval_f";
p.eval_u = "eval_u";
p.eval_fu = "F_U";
eval_linf = "eval_linf";
eval_y = "eval_y";

p.kappa = 5;
p.Nx = 10;
p.Ny = 10;
p.Nz = 4;
p.hx = 1;
p.hy = 1;
p.hz = 1;

p.magBx = 0;
p.magBy = 0;
% p.magBz = 0.5;
p.appliedBz = 0;
p.periodic_x = 0;
p.periodic_y = 0;
p.periodic_z = 0;

p = contruct_indices(p);

p.LPHIX = construct_LPHIXm(p);
p.LPHIY = construct_LPHIYm(p);
p.LPHIZ = construct_LPHIZm(p);

x = sparse(1:(p.Nx-1)*(p.Ny-1)*(p.Nz-1),1, 1);
y1 = sparse((p.Nx-1)*(p.Ny-1)*(p.Nz-1),1);
y2 = sparse((p.Nx-1)*(p.Ny-1)*(p.Nz-1),1);
y3 = sparse((p.Nx-1)*(p.Ny-1)*(p.Nz-1),1);

p.X0 = [x;y1;y2;y3];
p.U0 = 10;
p.magBz = p.U0;
p.t = 1;

p.t_start=0;
p.t_stop=10;
p.timestep = 1e-1;
visualize = 0;
p.visualizeSave = 0;
p.linearize = 1;
p.c = sparse(1:(p.Nx-1)*(p.Ny-1)*(p.Nz-1),1, 1/numel(x), numel(p.X0), 1);

[X,p] = Trapezoidal(p.eval_f,p.X0,p,p.eval_u,p.t_start,p.t_stop,p.timestep,visualize);



figure(1000)
hold on



%% linearized

p.U0 = 1;
p.magBz = p.U0;
p.t = 1;


%define F
[u,p] = feval(p.eval_u, p.t,p.X0,p);
F = @(X) feval(p.eval_f, X, p, u);

err = 1e-5;
fu = 0;
J0F = eval_num_jac(p.X0, p.eval_f, p, p.eval_u, p.t, err, fu);
p.A = J0F;

fu = 1;
J0U = eval_num_jac(p.U0, p.eval_fu, p, p.eval_u, 1, err, fu);
K0 = F(p.X0) - J0F*p.X0 - J0U*p.U0;
p.B = [K0 J0U];

p2.Breal = p.Breal;

p.c = sparse(1:(p.Nx-1)*(p.Ny-1)*(p.Nz-1),1, 1/numel(x), numel(p.X0), 1);

    
p.t_start=0;
p.t_stop=10;
p.timestep = 1e-1;
visualize = 0;
p.visualizeSave = 0;
% [X,p] = Trapezoidal_lin(eval_linf,p.X0,p,p.t_start,p.t_stop,p.timestep,visualize);

    
%% reduced
p2.A = p.A;
p2.B = p.B;
p2.c = p.c;

clear p
p.eval_f = "eval_f";
p.eval_u = "eval_u";
p.eval_fu = "F_U";

p.Breal = p2.Breal;


p.A = p2.A;
p.B = p2.B;
p.c = p2.c;



p.q = 4 * 3*3*2;

p.kappa = 5;
p.Nx = 4;
p.Ny = 4;
p.Nz = 3;
p.hx = 1;
p.hy = 1;
p.hz = 1;

p.magBx = 0;
p.magBy = 0;
% p.magBz = 0.5;
p.appliedBz = 0;
p.periodic_x = 0;
p.periodic_y = 0;
p.periodic_z = 0;

p = contruct_indices(p);

p.LPHIX = construct_LPHIXm(p);
p.LPHIY = construct_LPHIYm(p);
p.LPHIZ = construct_LPHIZm(p);

x = sparse(1:(p.Nx-1)*(p.Ny-1)*(p.Nz-1),1, 1);
y1 = sparse((p.Nx-1)*(p.Ny-1)*(p.Nz-1),1);
y2 = sparse((p.Nx-1)*(p.Ny-1)*(p.Nz-1),1);
y3 = sparse((p.Nx-1)*(p.Ny-1)*(p.Nz-1),1);

p.X0 = [x;y1;y2;y3];
p.U0 = 1;
p.t = 1;

   
    [V,D] = eigs(p.A, size(p.A,1));
    Bt = V'*p.B;
    ct = V'*p.c;

    q = p.q;
    D_diag = diag(D);
    [D_sort, i] = sort(D_diag);
    
   
    Aq = diag(D_diag(end-q+1:end));
    Bt_q = Bt(end-q+1:end,:);
    ct_q = ct(end-q+1:end);  
    
    p.A = Aq;
    p.B = Bt_q;
    p.c = sparse(1:(p.Nx-1)*(p.Ny-1)*(p.Nz-1),1, 1/numel(x), numel(p.X0), 1);
%     p.c = ct_q;

% 
% X = p.X0;
% U = p.U0;
% 
% F = feval(eval_linf, X, p, U);

% xeval = J0F\-F(X0);
% linf = A*X0 + B*[1;U0];
% f = F(xeval)
% imagesc(abs(f-linf))
% colorbar



%%
p.linearize = 1;
p.t_start=0;
p.t_stop=10;
p.timestep = 1e-1;
visualize = 0;
p.visualizeSave = 0;
[X,p] = Trapezoidal_lin(eval_linf,p.X0,p,p.t_start,p.t_stop,p.timestep,visualize);


%% 
% p.visualizeSave=1;
% visualizeNetworkX(X,p)


