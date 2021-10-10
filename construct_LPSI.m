function LPSI = construct_LPSI(y,Nx, Ny, Nz)
% can be used on and phi^{x,y,z} for u_{i,j,k}==psi_{i,j,k}
    LPSI = zeros(Nx+Ny+Nz);
    for k = 1:Nz
        for j = 1:Ny
            for i = 1:Nx
                m = Nx*(j-1)+Nx*Ny*(k-1);
                if i~=1 && i~=Nx
                    LPSI(i+m,i-1+m) = exp(-1i*y(i-1,j,k));
                    LPSI(i+m,i+m) = -2;
                    LPSI(i+m,i+1+m) = exp(1i*y(i,j,k));
                end
                if i==1
                    LPSI(i+m,i+m) = -1;
                    LPSI(i+m,i+1+m) = exp(1i*y(i,j,k));
                end
                if i==Nx
                    LPSI(i+m,i-1+m) = exp(-1i*y(i-1,j,k));
                    LPSI(i+m,i+m) = -1;
                end
            end
        end
    end
end