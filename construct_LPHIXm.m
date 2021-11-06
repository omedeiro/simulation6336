function LPHI = construct_LPHIXm(p)
hx = p.hx;
hy = p.hy;
hz = p.hz; 
kappa = p.kappa; 
Nx = p.Nx;
Ny = p.Ny;
Nz = p.Nz;

LPHI = sparse((Nx+1)*(Ny+1)*(Nz+1), (Nx+1)*(Ny+1)*(Nz+1));

mk = (p.Nx+1)*(p.Ny+1);
mj = (p.Nx+1);
m = p.M2; 


LPHI(sub2ind(size(LPHI),m,m-1)) =  - 2*(kappa^2/hy^2 + kappa^2/hz^2);


if Ny > 1
    LPHI(sub2ind(size(LPHI),m,m+mj)) = kappa^2/hy^2;
    LPHI(sub2ind(size(LPHI),m,m-mj)) = kappa^2/hy^2;
end

            
% k TERMS
if Nz > 1
    LPHI(sub2ind(size(LPHI),m,m+mk)) = kappa^2/hz^2;
    LPHI(sub2ind(size(LPHI),m,m-mk)) = kappa^2/hz^2;
end

end