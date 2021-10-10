function FPHIY = construct_FPHIY(x, y1, y2, y3, hx, hz, kappa, Nx, Ny, Nz)

    for k = 2 : Nz
        for j = 2 : Ny
            for i = 2 : Nx
                FPHIY_1(i,j,k) = (kappa^2/hz^2)*(-y3(i,j+1,k)+y3(i,j,k)+y3(i,j+1,k-1)-y3(i,j,k-1))...
                                +(kappa^2/hx^2)*(-y1(i,j+1,k)+y1(i,j,k)+y1(i-1,j+1,k)-y1(i-1,j,k))...
                                +imag(exp(-1i*y2(i,j,k))*conj(x(i,j,k))*x(i,j+1,k));
            
            end
        end
    end
    FPHIY = cube2column(FPHIY_1);
end                            

