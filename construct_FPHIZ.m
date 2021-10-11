function FPHIZ = construct_FPHIZ(x, y1, y2, y3, hx, hy, kappa, Nx, Ny, Nz)
    FPHIZ_1 = zeros(Nx+1,Ny+1,Nz+1);
    for k = 2 : Nz
        for j = 2 : Ny
            for i = 2 : Nx
                FPHIZ_1(i,j,k) =(kappa^2/hx^2)*(-y1(i,j,k+1)+y1(i,j,k)+y1(i-1,j,k+1)-y1(i-1,j,k))...
                                +(kappa^2/hy^2)*(-y2(i,j,k+1)+y2(i,j,k)+y2(i,j-1,k+1)-y2(i,j-1,k))...
                                +imag(exp(-1i*y3(i,j,k))*conj(x(i,j,k))*x(i,j,k+1));          
            
            end
        end
    end
    
    FPHIZ = cube2column(FPHIZ_1);
end                            

