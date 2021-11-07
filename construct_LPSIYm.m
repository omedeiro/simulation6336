function LPSI = construct_LPSIYm(y,p)
% can be used on and phi^{x,y,z} for u_{i,j,k}==psi_{i,j,k}LPSI = zeros(Nx*Ny*Nz);
%LPSI = zeros((Nx+1)*(Ny+1)*(Nz+1));
LPSI = sparse((p.Nx+1)*(p.Ny+1)*(p.Nz+1), (p.Nx+1)*(p.Ny+1)*(p.Nz+1));

% mk = (p.Nx+1)*(p.Ny+1);
mj = (p.Nx+1);
m = p.M2;
            
LPSI(p.L_m) = -2;
            
if p.Ny > 1
    LPSI(p.L_pmj) = exp(1i*y(m));
    LPSI(p.L_mmj) = exp(-1i*y(m-mj));
end


end