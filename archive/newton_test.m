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

p.LPHIX = construct_LPHIXm(p);
p.LPHIY = construct_LPHIYm(p);
p.LPHIZ = construct_LPHIZm(p);

x = sparse(1:(p.Nx-1)*(p.Ny-1)*(p.Nz-1),1, 1);
y1 = sparse((p.Nx-1)*(p.Ny-1)*(p.Nz-1),1);
y2 = sparse((p.Nx-1)*(p.Ny-1)*(p.Nz-1),1);
y3 = sparse((p.Nx-1)*(p.Ny-1)*(p.Nz-1),1);

x_start = [x;y1;y2;y3];
% 
% t = 0.1;
% 
% err_J = 1e-2;
% J = eval_num_jac(x_start, eval_f, p, eval_u, t, err_J)



u = feval(eval_u, 0);
errf = 1E-8;
errDeltax = 1E-8;
relDeltax = 1E-8;
MaxIter = 200;
visualize = true;
tolrGCR = 1E-8;
epsMF = 1E-3;

[X_field,converged,errf_k,errDeltax_k,relDeltax_k,iterations] = NewtonGCR(eval_f,x_start,p,u,errf, errDeltax, relDeltax, MaxIter, visualize, false, 0, tolrGCR, epsMF)


psi = X(1:(p.Nx-1)*(p.Ny-1)*(p.Nz-1));
        psi2 = column2cube(abs(psi).^2, (p.Nx-1), (p.Ny-1), (p.Nz-1));
        psi_surf = psi2(:,:,1);
        figure(2)
        s = surf(psi_surf)
        view(0,90)
        colorbar
        caxis([0 max(max(psi_surf))])

%% 
visualizeNetwork(1, X_field,p)
