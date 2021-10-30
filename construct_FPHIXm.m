function FPHIX = construct_FPHIXm(X, p)
colN = (p.Nx+1)*(p.Ny+1)*(p.Nz+1);
x = X(1:colN);
y1 = X(colN+1:2*colN);
y2 = X(2*colN+1:3*colN);
y3 = X(3*colN+1:4*colN);

FPHIX = sparse((p.Nx-1)*(p.Ny-1)*(p.Nz-1),1);

    for k = 2 : p.Nz
        for j = 2 : p.Ny
            for i = 2 : p.Nx
                mj = (p.Nx+1)*j-(p.Nx+1)*(j-1);
                mk = (p.Nx+1)*(p.Ny+1)*k-(p.Nx+1)*(p.Ny+1)*(k-1);
                m = i+(p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
                
                M = i+(p.Nx-1)*(j-1)+(p.Nx-1)*(p.Ny-1)*(k-1);

                FPHIX(M) = (p.kappa^2/p.hy^2)*(-y2(m+1)+y2(m)+y2(m+1-mj)-y2(m-mj))...
                                +(p.kappa^2/p.hz^2)*(-y3(m+1)+y3(m)+y3(m+1-mk)-y3(m-mk))...
                                +imag(exp(-1i*y1(m))*conj(x(m))*x(m+1));
            end
        end
    end
end                            
