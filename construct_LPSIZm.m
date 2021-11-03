function LPSI = construct_LPSIZm(y,p)
Nx = p.Nx;
Ny = p.Ny;
Nz = p.Nz;
LPSI = sparse((Nx+1)*(Ny+1)*(Nz+1), (Nx+1)*(Ny+1)*(Nz+1));

    mk = (p.Nx+1)*(p.Ny+1);
    mj = (p.Nx+1);
    m = p.M2;
    
            
            LPSI(sub2ind(size(LPSI),m,m)) = -2;

            if Nz > 1
                LPSI(sub2ind(size(LPSI),m,m+mk)) = exp(1i*y(m));
                LPSI(sub2ind(size(LPSI),m,m-mk)) = exp(-1i*y(m-mk));
            end

end