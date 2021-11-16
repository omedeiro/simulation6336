function FPHIX = construct_FPHIXm(x, y1, y2, y3, p)


mk = (p.Nx+1)*(p.Ny+1);
mj = (p.Nx+1);

m = p.M2; %internal in terms of N+1. 2:N
M = p.m2; %internal 1:N-1

FPHIX = sparse(M, 1, (p.kappa^2/p.hy^2).*(-y2(m+1)+y2(m)+y2(m+1-mj)-y2(m-mj))...
                    +(p.kappa^2/p.hz^2).*(-y3(m+1)+y3(m)+y3(m+1-mk)-y3(m-mk))...
                    +imag(exp(-1i*y1(m)).*conj(x(m)).*x(m+1)) );

end