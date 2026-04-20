function LPHI = construct_LPHIZm(p)
hx = p.hx;
hy = p.hy;
hz = p.hz; 
kappa = p.kappa; 

mk = (p.Nx+1)*(p.Ny+1);
mj = (p.Nx+1);
m = p.M2;

if p.Nz>1
    N_L = (p.Nx+1)*(p.Ny+1)*(p.Nz+1);
else
    N_L = (p.Nx+1)*(p.Ny+1);
end
LPHI = sparse(N_L, N_L);
        
            
% Diagonal
LPHI = LPHI + sparse(m, m, - 2*(kappa^2/hx^2 + kappa^2/hy^2), N_L, N_L);

% i TERMS
if p.Nx > 1
    LPHI = LPHI + sparse(m, m+1, kappa^2/hx^2, N_L, N_L);
    LPHI = LPHI + sparse(m, m-1, kappa^2/hx^2, N_L, N_L);
end

            
% k TERMS
if p.Ny > 1
    LPHI = LPHI + sparse(m, m+mj, kappa^2/hy^2, N_L, N_L);
    LPHI = LPHI + sparse(m, m-mj, kappa^2/hy^2, N_L, N_L);
end


end