function LPHI = construct_LPHIY(hx, hy, hz, kappa, Nx, Ny, Nz)
LPHI = zeros(Nx*Ny*Nz);
for k = 1:Nz
    for j = 1:Ny
        for i = 1:Nx
            mj = Nx*j-Nx*(j-1);
            mk = Nx*Ny*k-Nx*Ny*(k-1);
            m = Nx*(j-1)+Nx*Ny*(k-1);
            
            % Diagonal
            if i ~= 1 && i ~= Nx
                LPHI(i+m,i+m) = LPHI(i+m,i+m) - 2*(kappa^2/hz^2 + kappa^2/hx^2);
            end
            
            if i == 1 || i == Nx
                LPHI(i+m,i+m) = LPHI(i+m,i+m) -(kappa^2/hz^2 + kappa^2/hx^2);
            end
            
            % i TERMS
            if Nx > 1
                if i == 1 
                    LPHI(i+m,i+m+1) = LPHI(i+m,i+m+1) + kappa^2/hx^2;
                end     

                if i == Nx 
                    LPHI(i+m,i+m-1) = LPHI(i+m,i+m-1) + kappa^2/hx^2;
                end
                
                if i > 1 && i < Nx
                    LPHI(i+m,i+m+1) = LPHI(i+m,i+m+1) + kappa^2/hx^2;
                    LPHI(i+m,i+m-1) = LPHI(i+m,i+m-1) + kappa^2/hx^2;
                end

            end

            
            % k TERMS
            if Nz > 1
                if k == 1 
                    LPHI(i+m,i+m+mk) = LPHI(i+m,i+m+mk) + kappa^2/hz^2;
                end     

                if k == Nz
                    LPHI(i+m,i+m-mk) = LPHI(i+m,i+m-mk) + kappa^2/hz^2;
                end

                if k > 1 && k < Nz
                    LPHI(i+m,i+m+mk) = LPHI(i+m,i+m+mk) + kappa^2/hz^2;
                    LPHI(i+m,i+m-mk) = LPHI(i+m,i+m-mk) + kappa^2/hz^2;
                end
            end
       end
    end
end

end