function LPSI = construct_LPSIZm(y,p)
% Nx = p.Nx;
% Ny = p.Ny;
% Nz = p.Nz;
if p.Nz > 1
    N_L = (p.Nx+1)*(p.Ny+1)*(p.Nz+1);
else
    N_L = (p.Nx+1)*(p.Ny+1);
end
mk = (p.Nx+1)*(p.Ny+1);
% mj = (p.Nx+1);
m = p.M2;

LPSI = sparse(m, m, -2, N_L, N_L);

if p.Nz > 1
    LPSI = LPSI + sparse(m, m+mk, exp(1i*y(m)), N_L, N_L);
    LPSI = LPSI + sparse(m, m-mk, exp(-1i*y(m-mk)), N_L, N_L);
end

end