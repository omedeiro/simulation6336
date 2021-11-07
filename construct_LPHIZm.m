function LPHI = construct_LPHIZm(p)
hx = p.hx;
hy = p.hy;
hz = p.hz; 
kappa = p.kappa; 
Nx = p.Nx;
Ny = p.Ny;
Nz = p.Nz;
%LPHI = sparse((Nx+1)*(Ny+1)*(Nz+1), (Nx+1)*(Ny+1)*(Nz+1));

mk = (p.Nx+1)*(p.Ny+1);
mj = (p.Nx+1);
m = p.M2;
   
            
            
% Diagonal
LPHI(sub2ind(size(LPHI),m,m)) = - 2*(kappa^2/hx^2 + kappa^2/hy^2);

% i TERMS
if Nx > 1
    LPHI(sub2ind(size(LPHI),m,m+1)) = kappa^2/hx^2;
    LPHI(sub2ind(size(LPHI),m,m-1)) = kappa^2/hx^2;
end

            
% k TERMS
if Ny > 1
    LPHI(sub2ind(size(LPHI),m,m+mj)) = kappa^2/hy^2;
    LPHI(sub2ind(size(LPHI),m,m-mj)) = kappa^2/hy^2;
end


end