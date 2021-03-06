function LPSI = construct_LPSIXm(y,p)

if p.Nz > 1
    N_L = (p.Nx+1)*(p.Ny+1)*(p.Nz+1);
else
    N_L = (p.Nx+1)*(p.Ny+1);
end
m = p.M2;

LPSI = sparse(m, m-1, exp(-1i*y(m-1)), N_L, N_L);
LPSI = LPSI + sparse(m, m, -2, N_L, N_L);
LPSI = LPSI + sparse(m, m+1, exp(1i*y(m)), N_L, N_L);

end