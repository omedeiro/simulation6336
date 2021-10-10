function LPHI = construct_LPHIZ(hx, hy, hz, kappa, Nx, Ny, Nz)
LPHI = zeros(Nx*Ny*Nz);
for k = 2:Nz
    for j = 2:Ny
        for i = 2:Nx
            mj = Nx*j-Nx*(j-1);
            mk = Nx*Ny*k-Nx*Ny*(k-1);
            m = Nx*(j-1)+Nx*Ny*(k-1);
            
            % Diagonal
            if i ~= 1 && i ~= Nx
                LPHI(i+m,i+m) = LPHI(i+m,i+m) - 2*(kappa^2/hx^2 + kappa^2/hy^2);
            end
            
            if i == 2 || i == Nx
                LPHI(i+m,i+m) = LPHI(i+m,i+m) -(kappa^2/hx^2 + kappa^2/hy^2);
            end
            
            % i TERMS
            if Nx > 1
                if i == 2
                    LPHI(i+m,i+m+1) = LPHI(i+m,i+m+1) + kappa^2/hx^2;
                end     

                if i == Nx 
                    LPHI(i+m,i+m-1) = LPHI(i+m,i+m-1) + kappa^2/hx^2;
                end
                
                if i > 2 && i < Nx
                    LPHI(i+m,i+m+1) = LPHI(i+m,i+m+1) + kappa^2/hx^2;
                    LPHI(i+m,i+m-1) = LPHI(i+m,i+m-1) + kappa^2/hx^2;
                end

            end

            
            % k TERMS
            if Ny > 1
                if j == 2
                    LPHI(i+m,i+m+mj) = LPHI(i+m,i+m+mj) + kappa^2/hy^2;
                end     

                if j == Ny
                    LPHI(i+m,i+m-mj) = LPHI(i+m,i+m-mj) + kappa^2/hy^2;
                end

                if j > 2 && j < Ny
                    LPHI(i+m,i+m+mj) = LPHI(i+m,i+m+mj) + kappa^2/hy^2;
                    LPHI(i+m,i+m-mj) = LPHI(i+m,i+m-mj) + kappa^2/hy^2;
                end
            end
       end
    end
end

end