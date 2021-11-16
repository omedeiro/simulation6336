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

    Bx(m,t) = 1/(p.hy*p.hz) * (PHIY(p.M2B) - PHIY(p.M2B+mk) - PHIZ(p.M2B) + PHIZ(p.M2B+mj));

    By(m,t) = 1/(p.hz*p.hx) * (PHIZ(p.M2B) - PHIZ(p.M2B+1) - PHIX(p.M2B) + PHIX(p.M2B+mk));

    Bz(m,t) = 1/(p.hx*p.hy) * (PHIX(p.M2B) - PHIX(p.M2B+mj) - PHIY(p.M2B) + PHIY(p.M2B+1));

end

end
