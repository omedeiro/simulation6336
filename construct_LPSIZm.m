function LPSI = construct_LPSIZm(y,p)
Nx = p.Nx;
Ny = p.Ny;
Nz = p.Nz;
LPSI = sparse((Nx+1)*(Ny+1)*(Nz+1), (Nx+1)*(Ny+1)*(Nz+1));

for k = 2:Nz
    for j = 2:Ny
        for i = 2:Nx
            mj = (Nx+1)*j-(Nx+1)*(j-1);
            mk = (Nx+1)*(Ny+1)*k-(Nx+1)*(Ny+1)*(k-1);
            m = i+(Nx+1)*(j-1)+(Nx+1)*(Ny+1)*(k-1);
            
            LPSI(m,m) = -2;

            if Nz > 1
                LPSI(m,m+mk) = exp(1i*y(m));
                LPSI(m,m-mk) = exp(-1i*y(m-mk));
            end
        end
    end
end
end