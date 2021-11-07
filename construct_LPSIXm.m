function LPSI = construct_LPSIXm(y,p)
Nx = p.Nx;
Ny = p.Ny;
Nz = p.Nz;
%LPSI = sparse((Nx+1)*(Ny+1)*(Nz+1), (Nx+1)*(Ny+1)*(Nz+1));
mk = (p.Nx+1)*(p.Ny+1);
mj = (p.Nx+1);
m = p.M2;

LPSI(sub2ind(size(LPSI),m,m-1)) = exp(-1i*y(m-1));
LPSI(sub2ind(size(LPSI),m,m)) = -2;
LPSI(sub2ind(size(LPSI),m,m+1)) = exp(1i*y(m));
 
end