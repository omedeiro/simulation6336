function LPHI = construct_LPHIXm(p)
% hx = p.hx;
hy = p.hy;
hz = p.hz; 
kappa = p.kappa; 
% Nx = p.Nx;
% Ny = p.Ny;
% Nz = p.Nz;

N_L = (p.Nx+1)*(p.Ny+1)*(p.Nz+1);

mk = (p.Nx+1)*(p.Ny+1);
mj = (p.Nx+1);
m = p.M2;

LPHI = sparse(m, m, - 2*(kappa^2/hy^2 + kappa^2/hz^2), N_L, N_L);


if p.Ny > 1
    LPHI = LPHI + sparse(m, m+mj, kappa^2/hy^2, N_L, N_L);
    LPHI = LPHI + sparse(m, m-mj, kappa^2/hy^2, N_L, N_L);
end

            
% k TERMS
if p.Nz > 1
    LPHI = LPHI + sparse(m, m+mk, kappa^2/hz^2, N_L, N_L);
    LPHI = LPHI + sparse(m, m-mk, kappa^2/hz^2, N_L, N_L);
end

end