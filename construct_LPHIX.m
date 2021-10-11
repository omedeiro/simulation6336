function LPHI = construct_LPHIX(hx, hy, hz, kappa, Nx, Ny, Nz)
LPHI = zeros((Nx+1)*(Ny+1)*(Nz+1));
for k = 2:Nz
    for j = 2:Ny
        for i = 2:Nx
%             mj = Nx*j-Nx*(j-1);
%             mk = Nx*Ny*k-Nx*Ny*(k-1);
%             m = Nx*(j-1)+Nx*Ny*(k-1);

            mj = (Nx+1)*j-(Nx+1)*(j-1);
            mk = (Nx+1)*(Ny+1)*k-(Nx+1)*(Ny+1)*(k-1);
            m = (Nx+1)*(j-1)+(Nx+1)*(Ny+1)*(k-1);
            
            % Diagonal
%             if i ~= 2 && i ~= Nx
                LPHI(i+m,i+m) = LPHI(i+m,i+m) - 2*(kappa^2/hy^2 + kappa^2/hz^2);
%             end
            
%             if i == 2 || i == Nx
%                 LPHI(i+m,i+m) = LPHI(i+m,i+m) -(kappa^2/hy^2 + kappa^2/hz^2);
%             end
            
            % j TERMS
            if Ny > 1
%                 if j == 2 
%                     LPHI(i+m,i+m+mj) = LPHI(i+m,i+m+mj) + kappa^2/hy^2;
%                 end     
% 
%                 if j == Ny 
%                     LPHI(i+m,i+m-mj) = LPHI(i+m,i+m-mj) + kappa^2/hy^2;
%                 end
%                 
%                 if j > 2 && j < Ny
                    LPHI(i+m,i+m+mj) = LPHI(i+m,i+m+mj) + kappa^2/hy^2;
                    LPHI(i+m,i+m-mj) = LPHI(i+m,i+m-mj) + kappa^2/hy^2;
%                 end

            end

            
            % k TERMS
            if Nz > 1
%                 if k == 2
%                     LPHI(i+m,i+m+mk) = LPHI(i+m,i+m+mk) + kappa^2/hz^2;
%                 end     
% 
%                 if k == Nz
%                     LPHI(i+m,i+m-mk) = LPHI(i+m,i+m-mk) + kappa^2/hz^2;
%                 end

%                 if k > 2 && k < Nz
                    LPHI(i+m,i+m+mk) = LPHI(i+m,i+m+mk) + kappa^2/hz^2;
                    LPHI(i+m,i+m-mk) = LPHI(i+m,i+m-mk) + kappa^2/hz^2;
%                 end
            end
       end
    end
end

end