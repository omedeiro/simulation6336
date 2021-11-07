function LPSI = construct_LPSIZm(y,p)
% Nx = p.Nx;
% Ny = p.Ny;
% Nz = p.Nz;
LPSI = sparse((p.Nx+1)*(p.Ny+1)*(p.Nz+1), (p.Nx+1)*(p.Ny+1)*(p.Nz+1));

mk = (p.Nx+1)*(p.Ny+1);
% mj = (p.Nx+1);
m = p.M2;

            
LPSI(p.L_m) = -2;

if Nz > 1
    LPSI(p.L_pmk) = exp(1i*y(m));
    LPSI(p.L_mmk) = exp(-1i*y(m-mk));
end

end