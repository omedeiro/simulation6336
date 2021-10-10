function LPHI = construct_LPHIX(hx, hy, hz, kappa, Nx, Ny, Nz)
LPHI = zeros(Nx*Ny*Nz);
for k = 1:Nz
    for j = 1:Ny
        for i = 1:Nx
            m = Nx*(j-1)+Nx*Ny*(k-1);
            
            % Diagonal
            if i ~= 1 && i ~= Nx
                LPHI(i+m,i+m) = LPHI(i+m,i+m) - 2*(kappa^2/hy^2 + kappa^2/hz^2);
            end
            
            if i == 1 || i == Nx
                LPHI(i+m,i+m) = LPHI(i+m,i+m) -(kappa^2/hy^2 + kappa^2/hz^2);
            end
            
            
            if j ~= 1 && j ~= Ny && i~=1
                LPHI(i+m,i+m-1) = LPHI(i+m,i+m-1) + kappa^2/hy^2;
            end
            if j ~= 1 && j ~= Ny && i~=Nx
                LPHI(i+m,i+m+1) = LPHI(i+m,i+m+1) + kappa^2/hy^2;
            end
             
            if j == 1 && i ~= Nx
                LPHI(i+m,i+m+1) = LPHI(i+m,i+m+1) + kappa^2/hy^2;
            end     
             
            if j == 1 && i ~=1
                LPHI(i+m,i+m-1) = LPHI(i+m,i+m-1) + kappa^2/hy^2;
            end
            
            if j == Ny && i~=1
                LPHI(i+m,i+m-1) = LPHI(i+m,i+m-1) + kappa^2/hy^2;
            end
            
            if j == Ny && i ~= Nx
                LPHI(i+m,i+m+1) = LPHI(i+m,i+m+1) + kappa^2/hy^2;
            end   
            
       end
    end
end