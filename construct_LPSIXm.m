function LPSI = construct_LPSIXm(y,p)
% Nx = p.Nx;
% Ny = p.Ny;
% Nz = p.Nz;
LPSI = sparse((p.Nx+1)*(p.Ny+1)*(p.Nz+1), (p.Nx+1)*(p.Ny+1)*(p.Nz+1));
% mk = (p.Nx+1)*(p.Ny+1);
% mj = (p.Nx+1);

m = p.M2;

LPSI(p.L_m1) = exp(-1i*y(m-1));
LPSI(p.L_m) = -2;
LPSI(p.L_p1) = exp(1i*y(m));
 
end