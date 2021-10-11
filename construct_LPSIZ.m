function LPSI = construct_LPSIZ(y,Nx, Ny, Nz)
% can be used on and phi^{x,y,z} for u_{i,j,k}==psi_{i,j,k}LPSI = zeros(Nx*Ny*Nz);
LPSI = zeros((Nx+1)*(Ny+1)*(Nz+1));
for k = 2:Nz
    for j = 2:Ny
        for i = 2:Nx
%             mj = Nx*j-Nx*(j-1);
%             mk = Nx*Ny*k-Nx*Ny*(k-1);
%             m = Nx*(j-1)+Nx*Ny*(k-1);
            mj = (Nx+1)*j-(Nx+1)*(j-1);
            mk = (Nx+1)*(Ny+1)*k-(Nx+1)*(Ny+1)*(k-1);
            m = (Nx+1)*(j-1)+(Nx+1)*(Ny+1)*(k-1);
            
%             if i~=2 && i~=Nx
                LPSI(i+m,i+m) = -2;
%             end
%             if i==2 || i ==Nx
%                 LPSI(i+m,i+m) = -1;
%             end

            if Nz > 1
%                 if k == 2 
%                     LPSI(i+m,i+m+mk) = exp(1i*y(i,j,k));
%                 end     
% 
%                 if k == Nz 
%                     LPSI(i+m,i+m-mk) = exp(-1i*y(i,j,k-1));
%                 end
%                 
%                 if k > 2 && k < Nz
                    LPSI(i+m,i+m+mk) = exp(1i*y(i,j,k));
                    LPSI(i+m,i+m-mk) = exp(-1i*y(i,j,k-1));
%                 end

            end
        end
    end
end
end