function [Bx, By, Bz] = eval_Bfield(X,p)

mk = (p.Nx-1)*(p.Ny-1);
mj = (p.Nx-1);
% m = p.M2; 

n = (p.Nx-1)*(p.Ny-1)*(p.Nz-1);
n2 = (p.Nx-2)*(p.Ny-2)*(p.Nz-2);

Bx = sparse(n2,size(X,2));
By = sparse(n2,size(X,2));
Bz = sparse(n2,size(X,2));
m=1:n2;
for t = 1:size(X,2)
    PSI = X(1:n,t);
    PHIX = X(n+1:2*n,t);
    PHIY = X(2*n+1:3*n,t);
    PHIZ = X(3*n+1:4*n,t);

    Bx(m,t) = 1/(p.hy*p.hz) * (PHIY(m) - PHIY(m+mk) - PHIZ(m) + PHIZ(m+mj));

    By(m,t) = 1/(p.hz*p.hx) * (PHIZ(m) - PHIZ(m+1) - PHIX(m) + PHIX(m+mk));

    Bz(m,t) = 1/(p.hx*p.hy) * (PHIX(m) - PHIX(m+mj) - PHIY(m) + PHIY(m+1));


end
end
