function FPHIZ = construct_FPHIZ(x, y1, y2, y3, hx, hy, kappa, Nx, Ny, Nz)

    for k = 2 : Nz-1
        for j = 2 : Ny-1
            for i = 2 : Nx-1
                FPHIZ(i,j,k) =(kappa^2/hx^2)*(-y1(i,j,k+1)+y1(i,j,k)+y1(i-1,j,k+1)-y1(i-1,j,k))...
                                +(kappa^2/hy^2)*(-y2(i,j,k+1)+y2(i,j,k)+y2(i,j-1,k+1)-y2(i,j-1,k))...
                                +imag(exp(-1i*y3(i,j,k))*conj(x(i,j,k))*x(i,j,k+1));          
            
            end
        end
    end
end                            

