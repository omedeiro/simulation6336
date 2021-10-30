function LPHI = construct_LPHIZ(p)
hx = p.hx;
hy = p.hy;
hz = p.hz; 
kappa = p.kappa; 
Nx = p.Nx;
Ny = p.Ny;
Nz = p.Nz;
LPHI = zeros((Nx+1)*(Ny+1)*(Nz+1));
for k = 2:Nz
    for j = 2:Ny
        for i = 2:Nx
%             mj = Nx*j-Nx*(j-1);
%             mk = Nx*Ny*k-Nx*Ny*(k-1);
%             m = Nx*(j-1)+Nx*Ny*(k-1);
            mj = (Nx+1)*j-(Nx+1)*(j-1);
            mk = (Nx+1)*(Ny+1)*k-(Nx+1)*(Ny+1)*(k-1);
            m = i + (Nx+1)*(j-1)+(Nx+1)*(Ny+1)*(k-1);
            
            
            % Diagonal
%             if i ~= 1 && i ~= Nx
                LPHI(m,m) = LPHI(m,m) - 2*(kappa^2/hx^2 + kappa^2/hy^2);
%             end
            
%             if i == 2 || i == Nx
%                 LPHI(i+m,i+m) = LPHI(i+m,i+m) -(kappa^2/hx^2 + kappa^2/hy^2);
%             end
%             
            % i TERMS
            if Nx > 1
%                 if i == 2
%                     LPHI(i+m,i+m+1) = LPHI(i+m,i+m+1) + kappa^2/hx^2;
%                 end     
% 
%                 if i == Nx 
%                     LPHI(i+m,i+m-1) = LPHI(i+m,i+m-1) + kappa^2/hx^2;
%                 end
%                 
%                 if i > 2 && i < Nx
                    LPHI(m,m+1) = LPHI(m,m+1) + kappa^2/hx^2;
                    LPHI(m,m-1) = LPHI(m,m-1) + kappa^2/hx^2;
%                 end

            end

            
            % k TERMS
            if Ny > 1
%                 if j == 2
%                     LPHI(i+m,i+m+mj) = LPHI(i+m,i+m+mj) + kappa^2/hy^2;
%                 end     
% 
%                 if j == Ny
%                     LPHI(i+m,i+m-mj) = LPHI(i+m,i+m-mj) + kappa^2/hy^2;
%                 end

%                 if j > 2 && j < Ny
                    LPHI(m,m+mj) = LPHI(m,m+mj) + kappa^2/hy^2;
                    LPHI(m,m-mj) = LPHI(m,m-mj) + kappa^2/hy^2;
%                 end
            end
       end
    end
end

end