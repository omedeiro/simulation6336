function LPSI = construct_LPSIXm(y,p)
% Nx = p.Nx;
% Ny = p.Ny;
% Nz = p.Nz;
N_L = (p.Nx+1)*(p.Ny+1)*(p.Nz+1);
m = p.M2;

LPSI = sparse(m, m-1, exp(-1i*y(m-1)), N_L, N_L);
LPSI = LPSI + sparse(m, m, -2, N_L, N_L);
LPSI = LPSI + sparse(m, m+1, exp(1i*y(m)), N_L, N_L);

% LPSI(p.L_m1) = exp(-1i*y(m-1));
% LPSI(p.L_m) = -2;
% LPSI(p.L_p1) = exp(1i*y(m));
 
end