function LPHI = construct_LPHIZm(p)
hx = p.hx;
hy = p.hy;
hz = p.hz; 
kappa = p.kappa; 
Nx = p.Nx;
Ny = p.Ny;
Nz = p.Nz;
LPHI = sparse((Nx+1)*(Ny+1)*(Nz+1), (Nx+1)*(Ny+1)*(Nz+1));

    mk = (p.Nx+1)*(p.Ny+1);
    mj = (p.Nx+1);
    m = p.M2;
   
            
            
            % Diagonal
%             if i ~= 1 && i ~= Nx
                LPHI(sub2ind(size(LPHI),m,m)) = LPHI(sub2ind(size(LPHI),m,m)) - 2*(kappa^2/hx^2 + kappa^2/hy^2);
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
                    LPHI(sub2ind(size(LPHI),m,m+1)) = LPHI(sub2ind(size(LPHI),m,m+1)) + kappa^2/hx^2;
                    LPHI(sub2ind(size(LPHI),m,m-1)) = LPHI(sub2ind(size(LPHI),m,m-1)) + kappa^2/hx^2;
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
                    LPHI(sub2ind(size(LPHI),m,m+mj)) = LPHI(sub2ind(size(LPHI),m,m+mj)) + kappa^2/hy^2;
                    LPHI(sub2ind(size(LPHI),m,m-mj)) = LPHI(sub2ind(size(LPHI),m,m-mj)) + kappa^2/hy^2;
%                 end
            end


end