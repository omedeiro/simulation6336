function LPHI = construct_LPHIYm(p)
hx = p.hx;
% hy = p.hy;
hz = p.hz; 
kappa = p.kappa; 
% Nx = p.Nx;
% Ny = p.Ny;
% Nz = p.Nz;
mk = (p.Nx+1)*(p.Ny+1);
mj = (p.Nx+1);
m = p.M2;
N_L = (p.Nx+1)*(p.Ny+1)*(p.Nz+1);

            
% Diagonal

LPHI = sparse(m, m, - 2*(kappa^2/hz^2 + kappa^2/hx^2), N_L, N_L);

% i TERMS
if p.Nx > 1
    LPHI = LPHI + sparse(m, m+1, kappa^2/hx^2, N_L, N_L);
    LPHI = LPHI + sparse(m, m-1, kappa^2/hx^2, N_L, N_L);
end
   
% k TERMS
if p.Nz > 1
    LPHI = LPHI + sparse(m, m+mk, kappa^2/hz^2, N_L, N_L);
    LPHI = LPHI + sparse(m, m-mk, kappa^2/hz^2, N_L, N_L);
end


end