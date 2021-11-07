function LPHI = construct_LPHIXm(p)
% hx = p.hx;
hy = p.hy;
hz = p.hz; 
kappa = p.kappa; 
% Nx = p.Nx;
% Ny = p.Ny;
% Nz = p.Nz;

LPHI = sparse((p.Nx+1)*(p.Ny+1)*(p.Nz+1), (p.Nx+1)*(p.Ny+1)*(p.Nz+1));

% mk = (p.Nx+1)*(p.Ny+1);
% mj = (p.Nx+1);

LPHI(p.L_m) =  - 2*(kappa^2/hy^2 + kappa^2/hz^2);


if Ny > 1
    LPHI(p.L_pmj) = kappa^2/hy^2;
    LPHI(p.L_mmj) = kappa^2/hy^2;
end

            
% k TERMS
if Nz > 1
    LPHI(p.L_pmk) = kappa^2/hz^2;
    LPHI(p.L_mmk) = kappa^2/hz^2;
end

end