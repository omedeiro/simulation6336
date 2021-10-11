function LPSI = construct_LPSIX(y,Nx, Ny, Nz)
% can be used on and phi^{x,y,z} for u_{i,j,k}==psi_{i,j,k}
% Dimension of LPSIX should be ...
    LPSI = zeros((Nx+1)*(Ny+1)*(Nz+1));
    for k = 2:Nz
        for j = 2:Ny
            for i = 2:Nx
%                 m = Nx*(j-2)+Nx*Ny*(k-2);

%             mj = (Nx+1)*j-(Nx+1)*(j-1);
%             mk = (Nx+1)*(Ny+1)*k-(Nx+1)*(Ny+1)*(k-1);
                m = (Nx+1)*(j-1)+(Nx+1)*(Ny+1)*(k-1);
            
%                 m = m-1;
%                 if i~=2 && i~=Nx
                    LPSI(i+m,i-1+m) = exp(-1i*y(i-1,j,k));
                    LPSI(i+m,i+m) = -2;
                    LPSI(i+m,i+1+m) = exp(1i*y(i,j,k));
%                 end
%                 if i==2
%                     LPSI(i+m,i+m) = -1;
%                     LPSI(i+m,i+1+m) = exp(1i*y(i,j,k));
%                 end
%                 if i==Nx
%                     LPSI(i+m,i-1+m) = exp(-1i*y(i-1,j,k));
%                     LPSI(i+m,i+m) = -1;
%                 end
            end
        end
    end
end