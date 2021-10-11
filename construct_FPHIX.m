function FPHIX = construct_FPHIX(x, y1, y2, y3, hy, hz, kappa, Nx, Ny, Nz)
    FPHIX_1 = zeros(Nx-1,Ny-1,Nz-1);
    for k = 2 : Nz
        for j = 2 : Ny
            for i = 2 : Nx
                FPHIX_1(i-1,j-1,k-1) = (kappa^2/hy^2)*(-y2(i+1,j,k)+y2(i,j,k)+y2(i+1,j-1,k)-y2(i,j-1,k))...
                                +(kappa^2/hz^2)*(-y3(i+1,j,k)+y3(i,j,k)+y3(i+1,j,k-1)-y3(i,j,k-1))...
                                +imag(exp(-1i*y1(i,j,k))*conj(x(i,j,k))*x(i+1,j,k));
            end
        end
    end
    FPHIX = cube2column(FPHIX_1);
end                            
