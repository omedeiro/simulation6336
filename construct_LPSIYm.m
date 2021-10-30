function LPSI = construct_LPSIYm(y,p)
% can be used on and phi^{x,y,z} for u_{i,j,k}==psi_{i,j,k}LPSI = zeros(Nx*Ny*Nz);
% LPSI = zeros((Nx+1)*(Ny+1)*(Nz+1));
LPSI = sparse((p.Nx+1)*(p.Ny+1)*(p.Nz+1), (p.Nx+1)*(p.Ny+1)*(p.Nz+1));

for k = 2:p.Nz
    for j = 2:p.Ny
        for i = 2:p.Nx
            mj = (p.Nx+1)*j-(p.Nx+1)*(j-1);
            mk = (p.Nx+1)*(p.Ny+1)*k-(p.Nx+1)*(p.Ny+1)*(k-1);
            m = i+(p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
            
            LPSI(m,m) = -2;
            
            if p.Ny > 1
                LPSI(m,m+mj) = exp(1i*y(m));
                LPSI(m,m-mj) = exp(-1i*y(m-mj));
            end
        end
    end
end

end