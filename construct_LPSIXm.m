function LPSI = construct_LPSIX(y,p)
Nx = p.Nx;
Ny = p.Ny;
Nz = p.Nz;
LPSI = sparse((Nx+1)*(Ny+1)*(Nz+1), (Nx+1)*(Ny+1)*(Nz+1));
    for k = 2:Nz
        for j = 2:Ny
            for i = 2:Nx
                m = i+(Nx+1)*(j-1)+(Nx+1)*(Ny+1)*(k-1);
                LPSI(m,m-1) = exp(-1i*y(m-1));
                LPSI(m,m) = -2;
                LPSI(m,m+1) = exp(1i*y(m));
            end
        end
    end
end