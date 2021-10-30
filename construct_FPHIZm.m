function FPHIZ = construct_FPHIZm(X,p)
colN = (p.Nx+1)*(p.Ny+1)*(p.Nz+1);
x = X(1:colN);
y1 = X(colN+1:2*colN);
y2 = X(2*colN+1:3*colN);
y3 = X(3*colN+1:4*colN);

FPHIZ = sparse((p.Nx-1)*(p.Ny-1)*(p.Nz-1),1);
    for k = 2 : p.Nz
        for j = 2 : p.Ny
            for i = 2 : p.Nx
                mj = (p.Nx+1)*j-(p.Nx+1)*(j-1);
                mk = (p.Nx+1)*(p.Ny+1)*k-(p.Nx+1)*(p.Ny+1)*(k-1);
                m = i+(p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
                
                M = i+(p.Nx-1)*(j-1)+(p.Nx-1)*(p.Ny-1)*(k-1);

                
                FPHIZ(M) =(p.kappa^2/p.hx^2)*(-y1(m+mk)+y1(m)+y1(m-1+mk)-y1(m-1))...
                                +(p.kappa^2/p.hy^2)*(-y2(m+mk)+y2(m)+y2(m-mj+mk)-y2(m-mj))...
                                +imag(exp(-1i*y3(m))*conj(x(m))*x(m+mk));          
            
            end
        end
    end
end                            

