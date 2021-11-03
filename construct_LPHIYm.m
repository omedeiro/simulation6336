function LPHI = construct_LPHIYm(p)
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
 
%             mj = Nx*j-Nx*(j-1);
%             mk = Nx*Ny*k-Nx*Ny*(k-1);
%             m = Nx*(j-1)+Nx*Ny*(k-1);
            

            
            % Diagonal
%             if i ~= 2 && i ~= Nx
                LPHI(sub2ind(size(LPHI),m,m)) = LPHI(sub2ind(size(LPHI),m,m)) - 2*(kappa^2/hz^2 + kappa^2/hx^2);
%             end
            
%             if i == 2 || i == Nx
%                 LPHI(i+m,i+m) = LPHI(i+m,i+m) -(kappa^2/hz^2 + kappa^2/hx^2);
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
                
%                 if i > 2 && i < Nx
                    LPHI(sub2ind(size(LPHI),m,m+1)) = LPHI(sub2ind(size(LPHI),m,m+1)) + kappa^2/hx^2;
                    LPHI(sub2ind(size(LPHI),m,m-1)) = LPHI(sub2ind(size(LPHI),m,m-1)) + kappa^2/hx^2;
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
                    LPHI(sub2ind(size(LPHI),m,m+mk)) = LPHI(sub2ind(size(LPHI),m,m+mk)) + kappa^2/hz^2;
                    LPHI(sub2ind(size(LPHI),m,m-mk)) = LPHI(sub2ind(size(LPHI),m,m-mk)) + kappa^2/hz^2;
%                 end
            end


end